"""Software regression tests; fixtures do not validate dtangle or empirical FDR."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import next_core as N
import numpy as np
import pandas as pd
from common_support import intersection_budget
from matched_markers import allocate, freeze_matched, variants
from next_report import reliability, qualified_diagnostics, effective_counts, bh, validate_associations, FIELDS, CONTEXT
from clinical_probes import verify_baseline
from support import checkpoint


class BudgetTests(unittest.TestCase):
    def fixture(self):
        obs = np.ones((3, 14), int)
        obs[0, 3] = 0
        obs[1, 10] = 0
        owner = np.repeat(np.arange(1, 8), 2)
        weights = np.ones(14)
        sets = {'zero': [], 'legacy': np.arange(0, 14, 2), 'random_001': np.arange(1, 14, 2)}
        requested = np.full((3, 7), 8)
        return obs, owner, weights, sets, requested

    def test_same_mask_and_integer_budget(self):
        obs, owner, weights, sets, requested = self.fixture()
        budget, mask = intersection_budget(obs, owner, weights, sets, requested)
        self.assertEqual(int(budget.sum()), 152)
        self.assertEqual(budget[0, 1], 0)
        self.assertEqual(budget[1, 5], 0)
        for name, selected in sets.items():
            if name == 'zero':
                continue
            for rule in ('allow_zero', 'positive_only'):
                value, failures = allocate(weights, selected, owner, budget, obs, rule, 19)
                self.assertFalse(failures)
                actual = np.column_stack([value[:, owner == k].sum(1) for k in range(1, 8)])
                np.testing.assert_array_equal(actual, budget)
                if rule == 'positive_only':
                    self.assertEqual(int(value[obs == 0].sum()), 0)

    def test_all_zero_draw_is_retained(self):
        obs, owner, weights, sets, requested = self.fixture()
        obs[:] = 0
        budget, mask = intersection_budget(obs, owner, weights, sets, requested)
        self.assertFalse(budget.any())
        value, failures = allocate(weights, sets['legacy'], owner, budget, obs, 'positive_only', 5)
        self.assertFalse(failures)
        self.assertEqual(value.shape, obs.shape)
        self.assertFalse(value.any())

    def test_zero_weights_define_support(self):
        obs, owner, weights, sets, requested = self.fixture()
        weights[0] = 0
        with self.assertRaises(ValueError):
            intersection_budget(obs, owner, weights, sets, requested)

    def test_deterministic_and_no_input_changes(self):
        obs, owner, weights, sets, requested = self.fixture()
        copy = requested.copy()
        a, _ = intersection_budget(obs, owner, weights, sets, requested)
        b, _ = intersection_budget(obs, owner, weights, sets, requested)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(requested, copy)


class DiagnosticTests(unittest.TestCase):
    def scores(self):
        rows = []
        for i in range(20):
            for m in ('tmm', 'nnls_linear'):
                rows.append(dict(target_source='TEST', target_donor='d%02d'%i, method=m, outcome='null',
                    availability=1., score=(.3 if m == 'tmm' else (i+1)*1e-13), fdp=.2 if i % 2 else 0.,
                    status='FINITE_DIFFERENCE_SENSITIVITY_NONSMOOTH' if m == 'tmm' else 'LOCAL_SENSITIVITY',
                    stable_nuisance_rank=True,
                    derivative_agreement=.7 if m == 'tmm' else .001,
                    local_smoothness_warning=m == 'tmm'))
        return pd.DataFrame(rows)

    def test_nonsmooth_not_qualified(self):
        r = reliability(self.scores())
        self.assertFalse(r[r.method.eq('tmm')].derivative_qualified.any())
        self.assertTrue(r[r.method.eq('tmm')].score.isna().all())
        self.assertTrue(r[r.method.eq('tmm')].original_score.eq(.3).all())

    def test_numerical_noise_tied_in_auc(self):
        held, performance, counts = qualified_diagnostics(self.scores())
        self.assertTrue(held[held.method.eq('nnls_linear')].score.eq(0).all())
        finite = performance.auc.dropna()
        self.assertTrue(len(finite) > 0)
        self.assertTrue(finite.eq(.5).all())
        self.assertEqual(int(counts[counts.method.eq('tmm')].qualified_scores.iloc[0]), 0)

    def test_unavailable_baseline_relabel(self):
        frame = self.scores()
        frame.loc[1, 'availability'] = 0
        frame.loc[1, 'score'] = np.nan
        result = reliability(frame)
        self.assertEqual(result.loc[1, 'reliability'], 'UNAVAILABLE_BASE_FIT_OR_SCORE')


class ReportingTests(unittest.TestCase):
    def test_effective_donors_not_planned_donors(self):
        detail, mc = [], []
        for i in range(3):
            row = dict(unit_id='u'+str(i), target_donor='d'+str(i), base_subset='random_001',
                reference_source='R', target_source='T', subset_kind='random', support_rule='positive_only',
                outcome='null', method='tmm', planned_records=10, feasible_records=10 if i == 0 else 0)
            row.update({f: .1 if i == 0 else np.nan for f in FIELDS})
            detail.append(row)
            mc.append(dict({k: row[k] for k in ['unit_id', 'target_donor', 'base_subset', *CONTEXT]},
                           **{f+'__count': 10 if i == 0 else 0 for f in FIELDS}))
        result = effective_counts(pd.DataFrame(detail), pd.DataFrame(mc)).iloc[0]
        self.assertEqual(result.planned_donors, 3)
        self.assertEqual(result['fdp__contributing_donors'], 1)
        self.assertEqual(result['fdp__usable_repetition_subset_pairs'], 10)

    def test_bh_uses_full_family(self):
        np.testing.assert_allclose(bh([.001, .03, .2, 1]), [.004, .06, .26666666666666666, 1])

    def test_clinical_tables_require_all_keys(self):
        rows = []
        for ref in ('LIBD', 'SEAAD'):
            for method in ('tmm', 'tmm_exclude_factor', 'tmm_exclude_both', 'marker_denominator_only', 'marker_matched'):
                for gene, p, q in [('1', .01, .02), ('2', .7, .7)]:
                    rows.append(dict(reference=ref, method=method, gene_id=gene, p=p, q=q, available=True, beta=1., se=.3))
        frame = pd.DataFrame(rows)
        self.assertEqual(len(validate_associations(frame, ['1', '2'])), 10)
        frame.loc[0, 'q'] = .04
        with self.assertRaises(ValueError):
            validate_associations(frame, ['1', '2'])

    def test_baseline_drift_fails(self):
        old = pd.DataFrame(dict(gene_id=['1', '2'], available=[True, True], beta=[1., 2.], se=[.2, .3], p=[.1, .2], q=[.2, .2]))
        fit = {'beta': np.array([1., 2.]), 'se': np.array([.2, .3])}
        verify_baseline(['1', '2'], fit, np.ones(2, bool), np.array([.1, .2]), np.array([.2, .2]), old)
        fit['beta'][0] = 1.1
        with self.assertRaises(ValueError):
            verify_baseline(['1', '2'], fit, np.ones(2, bool), np.array([.1, .2]), np.array([.2, .2]), old)


class IntegrityTests(unittest.TestCase):
    def test_compact_review_has_self_contained_manifest(self):
        from finish import make_review
        from support import seal
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'run'
            stage = root/'clinical_multistep'
            (stage/'LIBD').mkdir(parents=True)
            N.save_json(stage/'scope.json', {'software_fixture': True})
            N.save_json(stage/'LIBD'/'not_exported.json', {'large_working_detail': True})
            cfg = dict(signature='fixture_review', units=[], settings=dict(rep_start=10000, repetitions=2))
            N.save_json(root/'protocol.json', cfg)
            seal(stage, ['scope.json', 'LIBD/not_exported.json'], cfg['signature'])
            make_review(root, cfg, 'SOFTWARE_FIXTURE_ONLY')
            N.verify_manifest(root/'review')
            self.assertFalse((root/'review'/'clinical_multistep'/'complete.json').exists())
            self.assertTrue((root/'review'/'clinical_multistep'/'server_receipt_identity.json').is_file())

    def test_output_overlap_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for output in (root, root/'new'):
                with self.assertRaises(ValueError):
                    N.assert_separate(output, root)
            N.assert_separate(root/'a', root/'b')

    def test_manifest_and_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'run'
            (root/'review').mkdir(parents=True)
            N.save_json(root/'protocol.json', dict(signature='0123456789fixture'))
            N.save_json(root/'review'/'x.json', {'software_fixture': True})
            N.write_manifest(root/'review')
            exported = N.export_review(root)
            self.assertEqual(exported.name, 'p3n_0123456789.zip')
            (root/'review'/'x.json').write_bytes(b'changed')
            with self.assertRaises(ValueError):
                N.export_review(root)

    def test_no_compact_export_as_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            N.save_json(root/'protocol.json', {'signature': 'wrong'})
            with self.assertRaises(ValueError):
                N.find_parent(root, root)


class PipelineFixture(unittest.TestCase):
    def test_clinical_multistep_and_missing_exports(self):
        import copy
        from types import SimpleNamespace
        import revision_driver as D
        import clinical
        import clinical_probes
        from next_report import clinical_exports
        from test_pipeline_fixture import bridge_double
        with tempfile.TemporaryDirectory(prefix='p3_next_clinical_') as tmp:
            root = Path(tmp)
            parent, output = root/'parent', root/'addon'
            parent.mkdir(); output.mkdir()
            first = D.synthetic_jobs(root/'prepared')[0]
            second = copy.deepcopy(first)
            first['source']['label'] = 'LIBD'; second['source']['label'] = 'SEAAD'
            n, g = 288, 16461
            rng = np.random.default_rng(119)
            x = np.r_[np.zeros(70), np.ones(218)]
            p = rng.dirichlet(np.ones(7)*5, n)
            counts = np.maximum(1, rng.poisson(50, size=(g, n)))
            counts[:700] = rng.poisson((p@first['base']['mean_counts'].T).T*30)
            counts = np.maximum(counts, 5)
            ids, samples = [str(i+1) for i in range(g)], ['s%d'%i for i in range(n)]
            raw = pd.DataFrame(counts, index=ids, columns=samples)
            meta = pd.DataFrame(dict(sample_id=samples, group=np.where(x, 'AD', 'Control'),
                age_numeric=rng.uniform(60, 90, n), sex=np.where(rng.random(n)>.5, 'F', 'M')))
            cov = np.column_stack([rng.normal(size=n), meta.sex.eq('M').to_numpy(float)])
            mapping = pd.DataFrame(dict(ensembl=first['source']['genes'], entrez=ids[:700], one_to_one=True))
            args = SimpleNamespace(output=parent, cache=root, offline=True, r_timeout=600)
            r = dict(executable='NOT_REAL_R', library='NOT_REAL_LIBRARY')
            original_cfg = dict(signature='SOFTWARE_FIXTURE_PARENT', settings=dict(derivative_step=.001, probe_steps=[.01, .05, .1]))
            cfg = dict(signature='SOFTWARE_FIXTURE_ADDON', settings=dict(
                clinical_derivative_steps=[.01, .001, .0001, .00001], clinical_probe_steps=[.01, .05, .1]))
            with patch('bulk_data.prepare', return_value=(raw, meta, mapping, dict(software_fixture=True))), \
                 patch('bulk_data.cohort', return_value=(meta, x, cov, meta)), patch.object(D, 'run_r', bridge_double):
                clinical.analyze(args, [first, second], r, original_cfg)
                (output/'corrected').mkdir()
                clinical_exports(parent, output/'corrected')
                self.assertEqual(len(pd.read_csv(output/'corrected'/'clinical_disease_associations.csv.gz')), 161610)
                args.output = output
                clinical_probes.run(args, parent, cfg, r)
                clinical_probes.run(args, parent, cfg, r)
            self.assertTrue(checkpoint(output/'clinical_multistep', cfg['signature']))
            records = pd.read_csv(output/'clinical_multistep'/'derivative_steps.csv')
            self.assertEqual(len(records), 40)
            self.assertEqual(set(records.derivative_step), {.01, .001, .0001, .00001})
            self.assertEqual(len(pd.read_csv(output/'clinical_multistep'/'finite_perturbations.csv')), 30)
            self.assertEqual(len(pd.read_csv(output/'clinical_multistep'/'baseline_identity.csv')), 10)

    def test_support_run_resume_and_summarize(self):
        import revision_driver as D
        import revision_core as C
        import common_support as S
        from selection import freeze
        from test_pipeline_fixture import bridge_double
        with tempfile.TemporaryDirectory(prefix='p3_next_fixture_') as tmp:
            root = Path(tmp)
            jobs = D.synthetic_jobs(root/'prepared')
            freeze(jobs, 1)
            freeze_matched(jobs, 1, root/'matching')
            units = C.plan_units(jobs, smoke=True)[:1]
            settings = dict(C.SEEDS, rep_start=10000, repetitions=2, subset_ids=variants(jobs[0]['stress_subsets']))
            cfg = dict(signature='SOFTWARE_FIXTURE_ONLY', settings=settings, units=units)
            N.save_json(root/'protocol.json', cfg)
            S.JOBS = jobs
            S.OPTIONS = dict(settings, output=root, signature=cfg['signature'],
                             r_info=dict(executable='NOT_R', library='NONE'), r_timeout=60)
            with patch.object(D, 'run_r', bridge_double):
                first = S.run_batch((units[0], 10000))
                second = S.run_batch((units[0], 10000))
                S.run_batch((units[0], 10001))
            self.assertFalse(first['resumed'])
            self.assertTrue(second['resumed'])
            self.assertEqual(first['rows'], 7*2*8)
            S.summarize(root, cfg)
            self.assertTrue(checkpoint(root/'common_support', cfg['signature']))
            summary = pd.read_csv(root/'common_support'/'method_summary.csv', keep_default_na=False, na_values=[''])
            self.assertEqual(set(summary.outcome), {'null', 'up025'})
            self.assertTrue(summary.design_feasibility.eq(1).all())
            self.assertTrue((root/'common_support'/'dose_retention.csv').is_file())
            batch = S.folder_for(root, (units[0], 10000))
            (batch/'metrics.csv.gz').write_bytes(b'corrupted')
            with self.assertRaises(ValueError):
                S.run_batch((units[0], 10000))


if __name__ == '__main__':
    unittest.main(verbosity=2)
