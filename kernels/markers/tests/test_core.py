import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from anchor_core import bh, design, hc3, metadata, pseudobulk
from driver import complete, validate_checkpoint


def sample_meta():
    return pd.DataFrame(dict(SampleID=['s'+str(i) for i in range(18)],
        Diagnosis=['AD']*11+['Control']*7,Age=[80,90,82,90,85,88]*3,
        Sex=['F','M']*9,Batch=[1,2,3]*6,PMI=[4.,5.,6.]*6,RIN=[7.,8.,9.]*6))


class Tests(unittest.TestCase):
    def test_design(self):
        d=sample_meta()
        n=pd.Series(100,index=d.SampleID)
        m,x,a,s=design(d,n,20,'primary')
        self.assertEqual(s['status'],'READY')
        self.assertEqual(s['n'],18)
        self.assertEqual(s['residual_df'],12)
        self.assertEqual(x.columns[1],'AD')

    def test_no_adaptive_cell_gate(self):
        d=sample_meta();n=pd.Series(100,index=d.SampleID)
        n.iloc[11:16]=19
        m,x,a,s=design(d,n,20,'primary')
        self.assertEqual(s['status'],'UNAVAILABLE')
        self.assertEqual(a.included.sum(),13)

    def test_covariate_missingness(self):
        d=sample_meta();d.loc[0,'PMI']=np.nan
        n=pd.Series(100,index=d.SampleID)
        self.assertEqual(len(design(d,n,20,'primary')[0]),18)
        self.assertEqual(len(design(d,n,20,'extended')[0]),17)

    def test_hc3_against_direct_sandwich(self):
        rng=np.random.default_rng(122)
        x=np.c_[np.ones(20),np.tile([0,1],10),rng.normal(size=20)]
        y=x@np.array([2.,.5,.3])+rng.normal(size=20)
        f=hc3(y,x)
        inv=np.linalg.inv(x.T@x);b=inv@x.T@y
        resid=y-x@b;h=np.diag(x@inv@x.T)
        v=inv@x.T@np.diag((resid/(1-h))**2)@x@inv
        self.assertAlmostEqual(f['beta'],b[1])
        self.assertAlmostEqual(f['se'],np.sqrt(v[1,1]))

    def test_rank_not_silently_repaired(self):
        x=np.c_[np.ones(12),np.repeat([0,1],6),np.repeat([0,1],6)]
        self.assertFalse(hc3(np.arange(12.),x)['available'])

    def test_bh_and_dose(self):
        np.testing.assert_allclose(bh([.01,.02,1]),[.03,.03,1])
        self.assertNotAlmostEqual(np.log2(1.5),.5)
        self.assertAlmostEqual(2**(-np.log2(1.5)),2/3)

    def test_checkpoint_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'a.csv').write_text('x\n1\n')
            complete(p,{'v':1},['a.csv'])
            self.assertTrue(validate_checkpoint(p,{'v':1}))
            (p/'a.csv').write_text('x\n2\n')
            self.assertFalse(validate_checkpoint(p,{'v':1}))

    def test_aggregation_matching_and_conservation(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'counts.h5'
            a=sparse.csc_matrix([[1,2,3,9],[4,5,6,10]])
            with h5py.File(path,'w') as f:
                g=f.create_group('matrix')
                for k,v in [('data',a.data),('indices',a.indices),('indptr',a.indptr),('shape',a.shape)]:g[k]=v
                g['barcodes']=np.array(['b2','b1','b3','unused'],dtype='S')
                g['features/id']=np.array(['ENSG1.2','ENSG2.3'],dtype='S')
                g['features/feature_type']=np.array(['Gene Expression']*2,dtype='S')
            cells=pd.DataFrame(dict(Barcode=['b1','b2','b3'],SampleID=['s1','s1','s2'],cell_type=['Astro']*3))
            genes,symbols,matrix,groups,audit=pseudobulk(path,cells)
            np.testing.assert_array_equal(matrix,[[3,3],[9,6]])
            self.assertEqual(audit['unused_unannotated_cells'],1)
            self.assertEqual(genes.tolist(),['ENSG1','ENSG2'])
            cells.loc[0,'Barcode']='absent'
            with self.assertRaisesRegex(ValueError,'missing'):pseudobulk(path,cells)

    def test_metadata_requires_donor_consistency(self):
        with tempfile.TemporaryDirectory() as d:
            m=sample_meta().iloc[:2].copy();m.Diagnosis=['AD','Control']
            m['Barcode']=['a','b'];m['Cell.Type']='ASC'
            path=Path(d)/'meta.csv';m.to_csv(path,index=False)
            self.assertEqual(len(metadata(path)[1]),2)
            m.loc[1,'SampleID']=m.loc[0,'SampleID'];m.to_csv(path,index=False)
            with self.assertRaisesRegex(ValueError,'Inconsistent'):metadata(path)


if __name__=='__main__':unittest.main()
