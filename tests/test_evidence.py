import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from pipeline import evidence, imaging_audit as imaging, probe_audit
from pipeline.io import read_json
from pipeline.publication import PUB, identity_check, load_frame


class EvidenceTests(unittest.TestCase):
    @staticmethod
    def reference_fixture():
        from scipy import sparse
        classes=['Astro','EndoMural','Excit','Inhib','Micro','OPC','Oligo']
        ids=['g%03d'%i for i in range(210)]
        meta=pd.DataFrame([(f'd{d}',c) for d in range(4) for c in classes],
                          columns=['donor_id','cell_type'])
        ix=np.arange(210)//30
        counts=1+100*(ix[:,None]==(np.arange(len(meta))%7))
        return sparse.csr_matrix(counts),ids,meta,classes

    def test_transfer_selection_uses_only_reference_and_available_ids(self):
        v,ids,meta,classes=self.reference_fixture()
        supported=ids[1:]
        a=imaging.transfer_panel(v,ids,meta,classes,supported,['target_only'])
        self.assertEqual(len(a['ids']),175)
        self.assertNotIn(ids[0],a['ids'])
        self.assertEqual(len(set(a['ids'])),len(a['ids']))
        self.assertTrue(0<a['anchor']<=1)
        # Different order of genes/libraries must not reassign any tied markers.
        b=imaging.transfer_panel(v[::-1,::-1],ids[::-1],meta.iloc[::-1],classes,supported,[])
        self.assertEqual(a['ids'],b['ids'])
        np.testing.assert_allclose(a['A'],b['A'])
        with self.assertRaises(ValueError): imaging.transfer_panel(v,ids,meta,classes,supported,['d0'])

    def test_imaging_distinguishes_original_panel_failure_from_transfer_fit(self):
        from scipy.io import mmwrite
        from pipeline.io import atomic_json
        v,ids,meta,classes=self.reference_fixture()
        libraries=pd.DataFrame(dict(library=['L1','L2'],block=['b1','b2'],donor=['T1','T2']))
        truth=pd.DataFrame([(b,c,.12,1) for b in libraries.block for c in imaging.CLASSES],
                           columns=['block','cell_type','measured','sections'])
        target=pd.DataFrame(np.ones((len(ids),2)),index=ids,columns=libraries.library)
        def r_fixture(cmd,check):
            folder=Path(cmd[4]); config=read_json(folder/'imaging_plan.json')
            self.assertEqual(len(config['count_rows']),175)
            pd.DataFrame([dict(sample_id=l,method=m,**dict.fromkeys(classes,1/7)) for l in libraries.library
                          for m in ['dtangle_marker','dtangle_tmm']]).to_csv(folder/'imaging_fractions.csv',index=False)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for ref in ('LIBD','SEAAD'):
                r=p/'cross/references'/ref; r.mkdir(parents=True)
                mmwrite(r/'reference.mtx',v); meta.to_csv(r/'cells.csv',index=False)
                pd.DataFrame(dict(gene_id=ids)).to_csv(r/'genes.csv',index=False)
                atomic_json(p/'clinical/GSE125583'/ref/'plan.json',dict(genes=ids,classes=classes,
                    marker_rows=list(range(1,8)),marker_owner=list(range(1,8))))
            with patch.object(imaging,'matched_inputs',return_value=(target,libraries,truth,{'expected_heldout_donors':['T1','T2']})),patch.object(imaging.subprocess,'run',side_effect=r_fixture):
                scores,_,audit,markers=imaging.compute(p,p/'clinical',p/'cross',p/'work','Rscript',[])
            self.assertEqual(set(scores.panel),{'reference_transfer'})
            self.assertEqual(set(markers.panel),{'original_clinical','reference_transfer'})
            self.assertEqual([x['status'] for x in audit['panels']],['insufficient_marker_support','estimable']*2)

    def test_ranked_probes_are_blind_to_sample_order_and_paired(self):
        raw=np.ones((31,4)); ids=['g%02d'%i for i in reversed(range(31))]
        raw[8:20,0]=40; raw[3,:]=0
        marks=[1,2]; tests=np.arange(3,32)
        p=probe_audit.design(raw,ids,tests,marks)
        self.assertEqual(p,probe_audit.design(raw[:,::-1],ids,tests,marks))
        self.assertEqual(len(p),16)
        for down,up in zip(p[1:-1:2],p[2:-1:2]):
            self.assertEqual(down['rows'],up['rows'])
            self.assertEqual(down['step'],-up['step'])
            self.assertFalse(set(marks)&set(up['rows']))
            if up['kind'].startswith('five'): self.assertEqual(len(up['rows']),5)
        tied=probe_audit.design(np.ones((31,4)),ids,tests,marks)
        self.assertEqual(tied[1]['rows'],[29])
        with self.assertRaises(ValueError): probe_audit.design(raw,ids,tests,[3])
        with self.assertRaises(ValueError): probe_audit.design(raw*0,ids,tests,marks)

    def test_identity_and_full_publication_design(self):
        check=identity_check(load_frame(PUB/'data','clinical_linear_checks.csv'))
        self.assertEqual((check['reference'],check['method'],check['log_count_step']),('SEAAD','tmm',.05))
        self.assertAlmostEqual(check['max_absolute_error'],5.176414852314792e-15,delta=1e-27)
        with self.assertRaises(ValueError): identity_check(pd.DataFrame({'exact_identity_max_error':[np.nan]}))
        s=load_frame(PUB/'data','TableS16.csv')
        self.assertEqual((len(s),s.method.nunique()),(96,4))
        self.assertTrue((s.datasets/s.donors==200).all())

    def test_donor_equal_errors_and_merged_classes(self):
        libs=pd.DataFrame({'library':['l1','l2','l3','l4'],'block':['b1','b1','b2','b3'],'donor':['d1','d1','d1','d2']})
        truth=pd.DataFrame([(b,c,v,1) for b,vs in [('b1',[.1,.2]),('b2',[.1,.2]),('b3',[.7,.3])]
                            for c,v in zip(['Astro','OligoOPC'],vs)],columns=['block','cell_type','measured','sections'])
        pred=pd.DataFrame([dict(sample_id=l,method=m,Astro=.1,EndoMural=.1,Excit=.2,Inhib=.2,Micro=.2,Oligo=.15,OPC=.05)
                           for l in libs.library for m in ['dtangle_marker','dtangle_tmm']])
        s,b=imaging.score(pred,libs,truth)
        scores=s[s.method=='dtangle_tmm'].set_index('cell_type')
        self.assertAlmostEqual(scores.at['Astro','mae'],.3)
        self.assertAlmostEqual(scores.at['OligoOPC','mae'],.05)
        self.assertEqual(scores.at['Astro','donors'],2)
        self.assertEqual(b.loc[b.block.eq('b1'),'libraries'].unique().tolist(),[2])
        with self.assertRaises(ValueError): imaging.score(pd.concat([pred,pred.iloc[[0]]]),libs,truth)
        with self.assertRaises(ValueError): imaging.score(pred.iloc[1:],libs,truth)

    def test_checkpoints_recover_and_reject_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            def failed(s,k):
                (k/'x').write_bytes(b'unfinished')
                raise RuntimeError('interruption')
            with self.assertRaises(RuntimeError): evidence.checkpoint(root,'a',[1],failed)
            self.assertFalse((root/'a').exists())
            d=evidence.checkpoint(root,'a',[1],lambda s,k:(k/'x').write_bytes(b'complete'))
            self.assertEqual(evidence.checkpoint(root,'a',[1],failed),d)
            with self.assertRaises(ValueError): evidence.checkpoint(root,'a',[2],failed)
            (d/'x').write_bytes(b'changed')
            with self.assertRaises(ValueError): evidence.checkpoint(root,'a',[1],failed)
            with self.assertRaises(ValueError): evidence.checkpoint(root,'b',[1],lambda s,k:None)

    def test_minimal_usable_example(self):
        with tempfile.TemporaryDirectory() as tmp:
            work=Path(tmp)/'sample'
            result=evidence.demo(work)
            self.assertEqual(result['family'],5)
            self.assertLess(result['equivalent_beta_max'],1e-10)
            self.assertGreater(result['aligned_beta_max'],.001)
            out=read_json(work/'result/audit.json')
            self.assertEqual(len(out['inputs']),3)
            with self.assertRaises(ValueError): evidence.demo(work)
