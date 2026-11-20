import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'frozen'/'engine'),str(HERE/'frozen')]
import numpy as np
import pandas as pd
from scipy import sparse
from matched_markers import allocate,match_indices,freeze_matched,stress,variants
from coupling_math import geometry_change,exact_change
from enhancement_experiment import nnls_prediction
from enhancement_report import method_differences,subset_comparisons,FIELDS
from report import validate_diagnostic


class RuntimeTests(unittest.TestCase):
    def test_same_binary_venv_still_requires_reentry(self):
        import runtime
        with tempfile.TemporaryDirectory() as name:
            base=Path(name)/'base_python'
            venv=Path(name)/'venv_python'
            base.write_bytes(b'fixture')
            venv.write_bytes(b'fixture')
            args=SimpleNamespace(python=venv,project_root=Path(name))
            with patch.object(runtime.sys,'executable',str(base)), \
                 patch.object(Path,'resolve',return_value=base), \
                 patch.dict(os.environ,{},clear=True), \
                 patch('runtime.subprocess.run',return_value=SimpleNamespace(returncode=0)), \
                 patch('runtime.subprocess.call',return_value=37) as call:
                self.assertEqual(runtime.choose_python(args,['--workers','36']),37)
                self.assertEqual(call.call_args.args[0][0],str(venv))
                self.assertEqual(call.call_args.kwargs['env']['P3_REVIEW_REEXEC'],'1')

    def test_valid_current_interpreter_does_not_reenter(self):
        import runtime
        args=SimpleNamespace(python=Path(sys.executable),project_root=HERE)
        with patch('runtime.subprocess.run',return_value=SimpleNamespace(returncode=0)), \
             patch('runtime.subprocess.call') as call:
            self.assertIsNone(runtime.choose_python(args,[]))
            call.assert_not_called()


class MatchingTests(unittest.TestCase):
    def fixture(self):
        owners=np.repeat(np.arange(1,8),10)
        features=np.tile(np.column_stack([np.arange(10)/10,np.arange(10)/20,np.ones(10)/2]),(7,1))
        selected=np.flatnonzero(np.tile([False,False,True,True,False,False,False,False,False,False],7))
        return owners,features,selected

    def test_matching_is_deterministic_and_class_preserving(self):
        owners,features,selected=self.fixture()
        genes=np.array(['g%03d'%i for i in range(len(owners))])
        a,rows=match_indices(features,owners,genes,selected,51)
        b,_=match_indices(features,owners,genes,selected,51)
        np.testing.assert_array_equal(a,b)
        self.assertFalse(set(a)&set(selected))
        self.assertEqual(len(rows),len(selected))
        for k in range(1,8):
            self.assertEqual(sum(owners[a]==k),sum(owners[selected]==k))

    def test_reference_changes_can_change_selection(self):
        owners,features,selected=self.fixture()
        genes=np.array(['g%d'%i for i in range(len(owners))])
        a,_=match_indices(features,owners,genes,selected,1)
        b,_=match_indices(features,owners,genes,selected,2)
        self.assertEqual(len(a),len(b))

    def test_budget_and_zero_support(self):
        owners,features,selected=self.fixture()
        observed=np.ones((4,len(owners)),dtype=int)
        budget=np.arange(28).reshape(4,7)
        weights=np.arange(1,len(owners)+1)
        a,missing=allocate(weights,selected,owners,budget,observed,'positive_only',4)
        self.assertFalse(missing)
        np.testing.assert_array_equal(np.column_stack([a[:,owners==k].sum(1) for k in range(1,8)]),budget)
        observed[1,selected[owners[selected]==1]]=0
        b,missing=allocate(weights,selected,owners,budget,observed,'positive_only',4)
        self.assertIsNone(b)
        self.assertEqual(missing,[(1,1,int(budget[1,0]))])
        c,missing=allocate(weights,selected,owners,budget,observed,'allow_zero',4)
        self.assertFalse(missing)
        self.assertEqual(c[1,owners==1].sum(),budget[1,0])

    def test_negative_counts_and_fractional_budget_rejected(self):
        owners,_,selected=self.fixture()
        with self.assertRaises(ValueError):
            allocate(np.ones(len(owners)),selected,owners,np.ones((2,7))*.5,
                     np.ones((2,len(owners))),'allow_zero',1)
        with self.assertRaises(ValueError):
            allocate(np.ones(len(owners)),selected,owners,np.ones((2,7)),
                     -np.ones((2,len(owners))),'allow_zero',1)

    def test_frozen_generator_integration(self):
        import revision_core as C
        import revision_driver as D
        from selection import freeze
        with tempfile.TemporaryDirectory() as tmp:
            jobs=D.synthetic_jobs(Path(tmp)/'prepared')
            freeze(jobs,1)
            manifest=freeze_matched(jobs,1,Path(tmp)/'matching')
            self.assertTrue(manifest)
            unit=C.plan_units(jobs,smoke=True)[0]
            job=jobs[unit['job_index']]
            base=C.draw_base(job,unit,0,C.SEEDS)
            hashes=[]
            for variant in variants(job['stress_subsets']):
                value,record=stress(job,base,variant)
                if variant!='zero':
                    hashes.append(record['budget_sha256'])
                if value is not None:
                    self.assertFalse(value[base['group']==0].any())
            self.assertEqual(len(set(hashes)),1)


class GeometryTests(unittest.TestCase):
    def test_same_space_and_nnls_scaling(self):
        rng=np.random.default_rng(3)
        A=rng.uniform(.01,1,(40,7))
        p=rng.dirichlet(np.ones(7),30)
        y=p@A.T
        a,_=nnls_prediction(A,y,np.ones(len(y)))
        b,_=nnls_prediction(A,y,np.exp(rng.normal(size=len(y))))
        np.testing.assert_allclose(a,b,atol=1e-10)
        x=rng.normal(size=len(p))
        geometry=geometry_change(x,p,p[:,[1,0,2,3,4,5,6]])
        self.assertLess(geometry['coefficient_weight_change'],1e-12)
        actual,expected=exact_change(x,rng.normal(size=(len(p),4)),p,p[:,[1,0,2,3,4,5,6]])
        np.testing.assert_allclose(actual,expected,atol=1e-12)

    def test_q_p_clinical_scope(self):
        from clinical import compare_fraction_scales
        rng=np.random.default_rng(81)
        q=rng.dirichlet(np.ones(7)*3,80)
        p=q/np.arange(1,8)
        p/=p.sum(1,keepdims=True)
        summary,table=compare_fraction_scales(rng.normal(size=(80,10)),rng.normal(size=80),
                                             rng.normal(size=(80,2)),q,p,np.arange(10))
        self.assertEqual(len(table),20)
        self.assertEqual(summary['planned_genes'],10)
        self.assertTrue(table.groupby('fraction_scale').size().eq(10).all())


class ReportTests(unittest.TestCase):
    def test_eta_separates_diagnostic_and_regular_work(self):
        from run_coupling import remaining_hours
        self.assertIsNone(remaining_hours({'diagnostic':1,'regular':9},
                                         {'diagnostic':[100,100,100],'regular':[]},36))
        self.assertAlmostEqual(remaining_hours({'diagnostic':1,'regular':9},
                                              {'diagnostic':[100,100,100],'regular':[10,10,10]},1),190/3600)

    def frame(self):
        rows=[]
        for ref in ['LIBD','SEAAD']:
            for donor in ['d1','d2','d3','d4','d5']:
                for rep in range(2):
                    for subset in ['legacy','random_001','matched_001']:
                        for support in ['allow_zero','positive_only']:
                            for method in ['tmm','marker_matched','oracle']:
                                rows.append(dict(unit_id=ref+donor,reference_source=ref,target_source='external',
                                    target_donor=donor,rep=rep,base_subset=subset,subset_kind=subset.split('_')[0],
                                    support_rule=support,subset_id=subset+'__'+support,outcome='null',method=method,
                                    intervention_feasible=not(rep==1 and subset=='matched_001'),
                                    **{f:(.2 if method=='tmm' else .1) for f in FIELDS}))
        frame=pd.DataFrame(rows)
        frame.loc[~frame.intervention_feasible,FIELDS]=np.nan
        return frame

    def test_unpaired_results_rejected(self):
        frame=self.frame()
        with self.assertRaises(ValueError):
            method_differences(frame.iloc[:-1])

    def test_infeasible_not_false_zero_fdr(self):
        frame=self.frame()
        diffs=method_differences(frame)
        self.assertTrue(diffs.loc[~diffs.common_design_feasible,'fdp_minus_tmm'].isna().all())
        paired=subset_comparisons(frame)
        self.assertTrue(paired.loc[~paired.common_design_feasible,'fdp_difference'].isna().all())

    def test_reference_duplicates_share_donor_fold(self):
        frame=self.frame()
        frame['score']=np.arange(len(frame))/1000
        heldout,performance=validate_diagnostic(frame)
        self.assertTrue(heldout.groupby(['target_source','target_donor']).fold.nunique().eq(1).all())
        self.assertFalse(performance.empty)


if __name__=='__main__':
    unittest.main()
