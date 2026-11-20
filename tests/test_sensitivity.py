import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.cli import parser
from pipeline.sensitivity import aligned, execute


class SensitivityTests(unittest.TestCase):
    def test_alignment_uses_identifiers(self):
        frame = pd.DataFrame(dict(sample_id=['b', 'a'], value=[2, 1]))
        self.assertEqual(aligned(frame, ['a', 'b'], 'test').value.tolist(), [1, 2])

    def test_missing_sample_rejected(self):
        with self.assertRaises(ValueError):
            aligned(pd.DataFrame(dict(sample_id=['a'])), ['a', 'b'], 'test')

    def test_duplicate_sample_rejected(self):
        with self.assertRaises(ValueError):
            aligned(pd.DataFrame(dict(sample_id=['a', 'a'])), ['a', 'b'], 'test')

    def fixture(self, root):
        rng = np.random.default_rng(157)
        n = 80
        ids = [f's{i}' for i in range(n)]
        x = rng.integers(0, 2, n)
        p = rng.dirichlet(np.ones(3), n)
        y = rng.normal(size=(n, 12)) + x[:, None] * .2
        pd.DataFrame(y, columns=[f'g{i}' for i in range(12)]).assign(sample_id=ids).to_csv(root / 'outcomes.csv', index=False)
        pd.DataFrame(dict(sample_id=ids, exposure=x, age=rng.normal(size=n))).to_csv(root / 'design.csv', index=False)
        fractions = pd.concat([pd.DataFrame(p, columns=['a', 'b', 'c']).assign(sample_id=ids, method=m) for m in ['base', 'same']])
        fractions.sample(frac=1, random_state=3).to_csv(root / 'fractions.csv', index=False)
        return parser().parse_args(['sensitivity', '--outcomes', str(root / 'outcomes.csv'),
            '--design', str(root / 'design.csv'), '--fractions', str(root / 'fractions.csv'),
            '--baseline', 'base', '--covariates', 'age', '--output', str(root / 'result')])

    def test_end_to_end_identity_and_full_family(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = self.fixture(root)
            execute(args)
            delta = pd.read_csv(root / 'result/design_comparison.csv').iloc[0]
            self.assertEqual(delta.beta_max, 0)
            self.assertEqual(delta.bh_status_changes, 0)
            self.assertEqual(delta.se_ratio_median, 1)
            self.assertEqual(len(pd.read_csv(root / 'result/gene_statistics.csv.gz')), 24)

    def test_nonestimable_exposure_keeps_family(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = self.fixture(root)
            design = pd.read_csv(args.design)
            design.exposure = 1
            design.to_csv(args.design, index=False)
            execute(args)
            result = pd.read_csv(root / 'result/gene_statistics.csv.gz')
            self.assertFalse(result.available.any())
            self.assertTrue(result.p.eq(1).all())
            self.assertTrue(result.q.eq(1).all())

    def test_invalid_fractions_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = self.fixture(root)
            frame = pd.read_csv(args.fractions)
            frame.loc[0, 'a'] = -1
            frame.to_csv(args.fractions, index=False)
            with self.assertRaises(ValueError):
                execute(args)

    def test_rounded_fractions_do_not_add_a_dimension(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args = self.fixture(root)
            frame = pd.read_csv(args.fractions)
            frame[['a', 'b', 'c']] = frame[['a', 'b', 'c']].round(7)
            frame.to_csv(args.fractions, index=False)
            execute(args)
            summary = pd.read_csv(root / 'result/method_summary.csv')
            self.assertTrue(summary.composition_rank.eq(2).all())


if __name__ == '__main__':
    unittest.main()
