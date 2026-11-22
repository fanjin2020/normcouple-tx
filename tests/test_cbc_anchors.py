"""Unit checks for fixed-anchor normalization and the prespecified decision."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from pipeline.cbc_anchors import anchor_denominators, decision, marker_prediction, prediction, tmm_denominator
from pipeline.cli import main
from pipeline.io import ROOT, read_json
from pipeline.runtime import check_enhancement_environment


class EnhancementEnvironmentTests(unittest.TestCase):
    pins = {'numpy': '2.5.3', 'pandas': '2.3.3', 'scipy': '1.18.1',
            'h5py': '3.16.0', 'matplotlib': '3.10.6'}

    def test_archived_requirements_reproduce_archived_specification(self):
        pins = dict(line.split('==', 1) for line in
                    (ROOT / 'requirements-archived.txt').read_text(encoding='utf-8').splitlines())
        spec = read_json(ROOT / 'publication/input_definitions/computational_specification.json')
        self.assertEqual(pins, spec['software'])

    def test_recorded_versions_pass(self):
        with patch('pipeline.runtime.platform.python_version', return_value='3.12.13'), \
             patch('pipeline.runtime.import_module', side_effect=lambda name: SimpleNamespace(__version__=self.pins[name])):
            self.assertEqual(check_enhancement_environment(), {'python': '3.12.13', **self.pins})

    def test_python_patch_level_is_checked(self):
        with patch('pipeline.runtime.platform.python_version', return_value='3.12.14'), \
             patch('pipeline.runtime.import_module', side_effect=lambda name: SimpleNamespace(__version__=self.pins[name])):
            with self.assertRaisesRegex(ValueError, r'python: expected 3\.12\.13, found 3\.12\.14'):
                check_enhancement_environment()

    def test_cbc_rejects_mismatched_numpy_before_analysis(self):
        versions = dict(self.pins, numpy='2.2.6')
        with TemporaryDirectory() as root, \
             patch('pipeline.runtime.platform.python_version', return_value='3.12.13'), \
             patch('pipeline.runtime.import_module', side_effect=lambda name: SimpleNamespace(__version__=versions[name])), \
             patch('pipeline.package.verify'), patch('pipeline.cbc_anchors.execute') as execute:
            output = Path(root) / 'new'
            with self.assertRaisesRegex(ValueError, r'numpy: expected 2\.5\.3, found 2\.2\.6'):
                main(['cbc-upgrade', '--pseudobulk', 'absent.npz', '--imaging-cache', 'libd',
                      '--bulk-cache', 'geo', '--output', str(output), '--allow-downloads'])
            execute.assert_not_called()
            self.assertFalse(output.exists())


class AnchorTests(unittest.TestCase):
    def test_tmm_matches_reference_implementation(self):
        a = np.array([[100,50,27,3,0,6,9,16,54,7],
                      [115,41,23,0,1,5,25,50,43,8],
                      [73,32,13,9,7,1,12,40,11,3]])
        d = np.array([400,510,370])
        expected = [1.1629336980854559, .9937680431319124, .8652866791642637]
        np.testing.assert_allclose(tmm_denominator(a,d)/d,expected,atol=2e-14,rtol=0)

    def test_marker_contrasts_match_dtangle(self):
        a = np.array([[20,4,8,3],[9,6,50,17]])
        r = np.array([[.07,.06,.001,.002],[.006,.003,.05,.045]])
        expected = [[.5762125633307871,.4237874366692130],[.1428870072196935,.8571129927803065]]
        np.testing.assert_allclose(marker_prediction(a,[800,1000],r,np.array([0,0,1,1])),
                                   expected,atol=2e-14,rtol=0)

    def test_weighted_median_and_unsupported_rows(self):
        count = [[3,6,12,50],[0,0,0,40]]
        np.testing.assert_equal(anchor_denominators(count,np.ones(4),[1,1,6,1]),[12,np.nan])
        np.testing.assert_equal(anchor_denominators(count,np.ones(4),[1,1,1,1]),[6,np.nan])

    def test_test_rows_cannot_enter_protected_methods(self):
        x=np.array([[50,99,12,23,11,22,100,20],[40,70,24,44,15,25,60,5]])
        model=dict(mark=np.array([0,1,2,3]),owners=np.array([0,0,1,1]),test=np.array([6]),
            anchor=np.array([4,5,7]),ref=np.ones((2,8))/8,abundance=np.array([.04,.09,.03]),weights=[2,1,1])
        a,da=prediction(x,x.sum(1),model)
        y=x.copy()
        y[1,6]+=130
        b,db=prediction(y,y.sum(1),model)
        for protected in ('complete','mrn','anchor_equal','anchor_weighted'):
            np.testing.assert_array_equal(da[protected],db[protected])
            np.testing.assert_array_equal(a[protected],b[protected])
        self.assertFalse(np.array_equal(da['tmm'],db['tmm']))

    def test_decision_is_a_joint_gate(self):
        summary=pd.DataFrame([dict(method='anchor_weighted',fdr_hi=.05,availability=1.)])
        contrasts=pd.DataFrame([dict(method='anchor_weighted',comparator='mrn',direction='Micro',
            probe='marker',outcome='sparse',fdr_hi=.01,coverage_lo=-.019,power_lo=.002)])
        self.assertTrue(decision(summary,contrasts,False)['eligible'])
        self.assertFalse(decision(summary,contrasts,True)['eligible'])
        contrasts.loc[0,'fdr_hi']=.021
        self.assertFalse(decision(summary,contrasts,False)['eligible'])


if __name__=='__main__':
    unittest.main()
