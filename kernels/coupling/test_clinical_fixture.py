"""Run the clinical adapter with generated inputs and real R; not clinical evidence."""
import copy
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'frozen'/'engine'),str(HERE/'frozen')]
import numpy as np
import pandas as pd
import revision_driver as D
import clinical
from runtime import run_r


class ClinicalAdapter(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('P3_TEST_RSCRIPT'),'Set P3_TEST_RSCRIPT and P3_TEST_R_LIBRARY')
    def test_generated_inputs_real_r(self):
        self.exercise_adapter(run_r,os.environ['P3_TEST_RSCRIPT'],os.environ['P3_TEST_R_LIBRARY'])

    def test_generated_inputs_software_bridge(self):
        from test_pipeline_fixture import bridge_double
        self.exercise_adapter(bridge_double,'NOT_REAL_R','NOT_REAL_LIBRARY')

    def exercise_adapter(self,executor,executable,library):
        with tempfile.TemporaryDirectory(prefix='p3_clinical_test_') as name:
            root=Path(name)
            first=D.synthetic_jobs(root/'prepared')[0]
            second=copy.deepcopy(first)
            first['source']['label']='LIBD'; second['source']['label']='SEAAD'
            n,g=288,16461
            rng=np.random.default_rng(119)
            x=np.r_[np.zeros(70),np.ones(218)]
            p=rng.dirichlet(np.ones(7)*5,n)
            mu=first['base']['mean_counts']
            counts=np.maximum(1,rng.poisson(50,size=(g,n)))
            counts[:700]=rng.poisson((p@mu.T).T*30)
            counts=np.maximum(counts,5)
            ids=[str(i+1) for i in range(g)]
            samples=['s%d'%i for i in range(n)]
            raw=pd.DataFrame(counts,index=ids,columns=samples)
            meta=pd.DataFrame(dict(sample_id=samples,group=np.where(x,'AD','Control'),
                age_numeric=rng.uniform(60,90,n),sex=np.where(rng.random(n)>.5,'F','M')))
            cov=np.column_stack([rng.normal(size=n),meta.sex.eq('M').to_numpy(float)])
            mapping=pd.DataFrame(dict(ensembl=first['source']['genes'],entrez=ids[:700],one_to_one=True))
            args=SimpleNamespace(output=root,cache=root,offline=True,r_timeout=600)
            r=dict(executable=executable,library=library)
            cfg=dict(signature='fixture-only',settings=dict(derivative_step=.001,probe_steps=[.01,.05,.10]))
            with patch('bulk_data.prepare',return_value=(raw,meta,mapping,dict(software_fixture=True))), \
                 patch('bulk_data.cohort',return_value=(meta,x,cov,meta)),patch.object(D,'run_r',executor), \
                 patch.dict(os.environ,{'R_LIBS':library}):
                clinical.analyze(args,[first,second],r,cfg)
            stable=pd.read_csv(root/'clinical'/'method_stability.csv')
            self.assertEqual(len(stable),10)
            self.assertTrue(stable.planned_genes.eq(16161).all())
            inv=pd.read_csv(root/'clinical'/'invariance.csv')
            self.assertTrue(inv.max_change.lt(1e-10).all())
            assoc=pd.read_csv(root/'clinical'/'disease_associations.csv.gz')
            self.assertEqual(len(assoc),161610)
            self.assertTrue(assoc.available.all())
            scale=pd.read_csv(root/'clinical'/'clinical_fraction_scale.csv')
            self.assertEqual(len(scale),10)
            self.assertTrue(scale.planned_genes.eq(16161).all())
            detail=pd.read_csv(root/'clinical'/'fraction_scale_associations.csv.gz')
            self.assertEqual(len(detail),2*161610)
            self.assertEqual(set(detail.fraction_scale),{'rna_contribution','nuclear_fraction'})
            q=pd.read_csv(root/'clinical'/'rna_predictions.csv')
            self.assertEqual(len(q),288*10)
            self.assertEqual(set(q.fraction_definition),{'marker_RNA_contribution','full_reference_RNA_contribution'})


if __name__=='__main__':
    unittest.main()
