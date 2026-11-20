import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from pipeline.cli import parser
from pipeline.cross_method import paired_metrics, paired_precision, METHODS, PROBES, _valid_batch, _reference
from pipeline.io import atomic_json, receipt


class CrossMethodTests(unittest.TestCase):
    def test_cli_default_is_bounded(self):
        args = parser().parse_args(['cross-method', '--inputs', 'inputs', '--work', 'work', '--r-library', 'lib'])
        self.assertEqual(args.workers, 4)
        self.assertEqual(args.repetitions, 200)
        self.assertFalse(args.resume)

    def test_precision_uses_paired_differences(self):
        rows = []
        for method in METHODS:
            for donor in ('a','b'):
                for rep in range(4):
                    for probe in PROBES:
                        value = rep*.1 + (0 if probe=='baseline' else .02)
                        rows.append(dict(reference_source='R',target_source='T',target_donor=donor,
                            rep=rep,effect=0,method=method,probe=probe,fdp=value,power=np.nan,
                            coverage=value,proportion_rmse=value))
        summary = paired_precision(pd.DataFrame(rows))
        fdr = summary[summary.metric.eq('fdp')]
        self.assertTrue(np.allclose(fdr.donor_equal_change,.02))
        self.assertTrue(np.all(fdr.conditional_mc_se<1e-15))
        self.assertTrue(fdr.target_met.all())
        self.assertTrue(summary[summary.metric.eq('power')].donor_equal_change.isna().all())

    def test_incomplete_precision_pairs_rejected(self):
        rows = []
        for method in METHODS:
            for rep in range(2):
                for probe in PROBES:
                    rows.append(dict(reference_source='R',target_source='T',target_donor='d',rep=rep,
                        effect=0,method=method,probe=probe,fdp=0.,power=0.,coverage=1.,proportion_rmse=.1))
        with self.assertRaisesRegex(ValueError,'Incomplete paired'):
            paired_precision(pd.DataFrame(rows).iloc[:-1])

    def test_fixed_response_identity(self):
        rng = np.random.default_rng(41)
        y = rng.normal(size=(24, 8))
        w0, w1 = rng.normal(size=(2, 24))
        fit = lambda w: dict(weights=w, beta=w @ y, se=np.ones(8))
        props = np.full((24, 7), 1 / 7)
        result = paired_metrics(fit(w0), fit(w1), y, props, props)
        self.assertLess(result['weight_identity_error'], 1e-12)
        self.assertGreater(result['weight_l2'], 0)

    def test_invariant_prediction(self):
        rng = np.random.default_rng(42)
        y = rng.normal(size=(24, 8))
        w = rng.normal(size=24)
        fit = dict(weights=w, beta=w @ y, se=np.ones(8))
        props = np.full((24, 7), 1 / 7)
        result = paired_metrics(fit, fit, y, props, props)
        self.assertEqual(result['weight_l2'], 0)
        self.assertEqual(result['beta_rms'], 0)
        self.assertEqual(result['se_ratio_median'], 1)

    def test_missing_receipt_is_not_complete(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            atomic_json(root / 'plan.json', {})
            self.assertFalse(_valid_batch(root, 'abc'))

    def test_changed_batch_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            atomic_json(root / 'output.json', {'value': 1})
            atomic_json(root / 'receipt.json', {'protocol': 'abc', 'files': receipt(root, ['output.json'])})
            self.assertTrue(_valid_batch(root, 'abc'))
            atomic_json(root / 'output.json', {'value': 2})
            with self.assertRaises(ValueError):
                _valid_batch(root, 'abc')

    def test_mismatched_protocol_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            atomic_json(root / 'receipt.json', {'protocol': 'old', 'files': {'file': 'hash'}})
            with self.assertRaises(ValueError):
                _valid_batch(root, 'new')

    def test_interrupted_reference_is_rebuilt_then_sealed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder = root / 'references/LIBD'
            folder.mkdir(parents=True)
            (folder / 'reference.mtx').write_text('interrupted')
            source = dict(label='LIBD',x=sparse.csc_matrix([[1,2],[3,4]]),genes=np.array(['a','b']),
                meta=pd.DataFrame(dict(donor_id=['d1','d2'],cell_type=['A','A'])),
                masks={'d1':[np.array([0])],'d2':[np.array([1])]})
            job = dict(source=source, refs=['d1','d2'])
            _reference(job,root,{'LIBD/counts.mtx':'pinned'})
            self.assertTrue((folder/'source_receipt.json').is_file())
            _reference(job,root,{'LIBD/counts.mtx':'pinned'})
            with self.assertRaisesRegex(ValueError,'Reference inputs'):
                _reference(job,root,{'LIBD/counts.mtx':'different'})
            (folder/'cells.csv').write_text('changed')
            with self.assertRaises(ValueError):
                _reference(job,root,{'LIBD/counts.mtx':'pinned'})


if __name__ == '__main__':
    unittest.main()
