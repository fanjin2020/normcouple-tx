import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from pipeline.cli import parser
from pipeline.cohorts import expression_family, external_cohort, gene_mapping, soft_samples
from pipeline.clinical_estimators import METHODS, PROBES, summarize_fits, reuse_estimates
from pipeline.cross_method import _setup_runtime
from pipeline.run import commands, framework_identity, validate_run
from pipeline.io import atomic_json,digest,receipt,read_json


class ClinicalTests(unittest.TestCase):
    def test_mapping_is_one_to_one(self):
        frame=gene_mapping()
        self.assertTrue(frame.index.is_unique and frame.ensembl.is_unique)

    def test_family_excludes_all_markers_without_filtering_by_disease(self):
        raw=pd.DataFrame(np.ones((140,6)),index=[str(i) for i in range(140)])
        raw.iloc[0]=0
        selected=expression_family(raw,{'1','2'})
        self.assertEqual(selected.tolist(),list(range(3,140)))

    def test_external_uses_raw_fractional_columns_not_rpkm(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            text=[]
            data={'GeneID':[str(i) for i in range(101,221)],'GeneSymbol':['G']*120}
            for i in range(1,18):
                code=('C' if i<=8 else 'A')+str(i)
                status='control' if i<=8 else "advanced Alzheimer's Disease"
                text += ['^SAMPLE = GSM'+str(i),'!Sample_description = '+code,
                    '!Sample_characteristics_ch1 = disease status: '+status,
                    '!Sample_characteristics_ch1 = tissue subtype: dorsolateral prefrontal cortex']
                data[code+'_raw']=[i+.125]*120
                data[code+'_rpkm']=[999.]*120
            with gzip.open(root/'GSE53697_family.soft.gz','wt') as stream:
                stream.write('\n'.join(text))
            pd.DataFrame(data).to_csv(root/'GSE53697_RNAseq_AD.txt.gz',sep='\t',index=False)
            with patch('pipeline.cohorts.external_sources',return_value={}):
                raw,meta,provenance=external_cohort(root)
            self.assertEqual(raw.shape,(120,17))
            self.assertEqual(raw.iloc[0,0],1.125)
            self.assertEqual(provenance['fractional_entries'],120*17)
            self.assertEqual(meta.group.eq('AD').sum(),9)
            self.assertEqual(meta.sequencing_batch.eq(2).sum(),6)
            self.assertEqual(provenance['covariates'],['sequencing_batch'])

    def test_duplicate_metadata_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'x.gz'
            with gzip.open(path,'wt') as stream:
                stream.write('^SAMPLE = GSM1\n!Sample_characteristics_ch1 = sex: M\n!Sample_characteristics_ch1 = sex: F\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                soft_samples(path)

    def test_fixed_clinical_outcomes_produce_zero_change_for_identical_fractions(self):
        with tempfile.TemporaryDirectory() as temp:
            c,_,s=_setup_runtime(Path(temp))
            rng=np.random.default_rng(312)
            p=rng.dirichlet(np.ones(7),size=40)
            x=np.tile([0.,1.],20)
            y=rng.normal(size=(40,120))
            ids=['s'+str(i) for i in range(40)]
            frames=[]
            for method in METHODS:
                for probe in PROBES:
                    frame=pd.DataFrame(p,columns=c.M.CT)
                    frames.append(frame.assign(sample_id=ids,method=method,probe=probe))
            sm,cm,detail=summarize_fits(s,c,y,x,None,pd.concat(frames),c.M.CT,list(range(120)),ids)
            self.assertEqual(len(sm),18)
            self.assertEqual(len(cm),17)
            self.assertTrue(cm.bh_crossings.eq(0).all())
            self.assertTrue(cm.beta_max.eq(0).all())
            self.assertEqual(len(detail),18*120)

    def test_acceptance_does_not_export_as_research(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            cfg=dict(framework_source=framework_identity(),profile='acceptance')
            atomic_json(root/'pipeline.json',dict(cfg,signature=digest(cfg)))
            with self.assertRaisesRegex(ValueError,'not a full scientific'):
                validate_run(root)

    def test_unsupported_markers_do_not_block_native_music_or_count_as_zero_discoveries(self):
        with tempfile.TemporaryDirectory() as temp:
            c,_,s=_setup_runtime(Path(temp))
            rng=np.random.default_rng(711)
            p=rng.dirichlet(np.ones(7),size=40)
            x=np.tile([0.,1.],20)
            y=rng.normal(size=(40,120))
            ids=['s'+str(i) for i in range(40)]
            frames=[]
            for method in ('music_broad','music_non_test'):
                for probe in PROBES:
                    frames.append(pd.DataFrame(p,columns=c.M.CT).assign(sample_id=ids,method=method,probe=probe))
            predictions=pd.concat(frames)
            missing={m:'INSUFFICIENT_MARKER_COVERAGE' for m in METHODS if m not in ('music_broad','music_non_test')}
            with self.assertRaisesRegex(ValueError,'Missing clinical'):
                summarize_fits(s,c,y,x,None,predictions,c.M.CT,list(range(120)),ids)
            sm,cm,detail=summarize_fits(s,c,y,x,None,predictions,c.M.CT,list(range(120)),ids,missing)
            self.assertEqual(len(sm),18)
            self.assertEqual(len(cm),17)
            self.assertEqual(cm.available.sum(),4)
            unsupported=sm[sm.method.isin(missing)]
            self.assertTrue(unsupported.discoveries.isna().all())
            self.assertTrue(unsupported.available_tests.eq(0).all())
            self.assertEqual(len(detail),6*120)

    def test_single_pipeline_contains_new_clinical_stage(self):
        args=parser().parse_args(['run','--inputs','inputs','--work','work','--r-library','R',
                                 '--external-bulk-cache','external','--profile','acceptance'])
        plan=dict(commands(args,Path('runtime')))
        self.assertIn('clinical_estimators',plan)
        self.assertIn('--smoke',plan['cross_method'])
        self.assertIn('--smoke',plan['common_support'])
        self.assertNotIn('--synthetic',sum(plan.values(),[]))

    def test_unit_leverage_is_missing_inference_not_zero_discovery(self):
        with tempfile.TemporaryDirectory() as temp:
            c,_,s=_setup_runtime(Path(temp))
            rng=np.random.default_rng(917)
            p=np.zeros((40,7))
            p[1:,1:]=rng.dirichlet(np.ones(6),size=39)
            p[0,0]=1
            x=np.tile([0.,1.],20)
            y=rng.normal(size=(40,120))
            ids=['s'+str(i) for i in range(40)]
            frames=[pd.DataFrame(p,columns=c.M.CT).assign(sample_id=ids,method=m,probe=z)
                    for m in METHODS for z in PROBES]
            sm,cm,_=summarize_fits(s,c,y,x,None,pd.concat(frames),c.M.CT,list(range(120)),ids)
            self.assertTrue(sm.available_tests.eq(0).all())
            self.assertTrue(sm.discoveries.isna().all())
            self.assertTrue(sm.reason.eq('HC3_not_estimable').all())
            self.assertTrue(cm.bh_crossings.isna().all())
            self.assertFalse(cm.available.any())

    def test_estimator_cache_requires_sealed_outputs_identical_inputs_and_plan(self):
        with tempfile.TemporaryDirectory() as temp:
            source,dest=Path(temp)/'source',Path(temp)/'dest'
            old=source/'GSE53697/LIBD'
            old.mkdir(parents=True)
            dest.mkdir()
            prior=dict(schema=2,methods=list(METHODS),code={
                'pipeline/clinical_estimators.py':'old-report',
                'pipeline/clinical_estimators.R':'same-estimator'})
            current=dict(prior,code=dict(prior['code'],**{'pipeline/clinical_estimators.py':'new-report'}))
            atomic_json(source/'protocol.json',prior)
            for folder in (old,dest):
                (folder/'raw.bin').write_bytes(b'raw')
                (folder/'A.bin').write_bytes(b'reference')
                atomic_json(folder/'plan.json',dict(reference=str(folder/'bank'),genes=['g1'],test_rows=[1]))
            names=['fractions.csv','engine.json','outcomes.bin','outcome_denominator.csv','plan.json']
            for name in names[:-1]:
                (old/name).write_bytes(b'verified output')
            atomic_json(old/'receipt.json',dict(protocol=digest(prior),files=receipt(old,names)))
            atomic_json(source/'complete.json',dict(protocol=digest(prior),files=receipt(source,['GSE53697/LIBD/receipt.json'])))
            reuse_estimates(source,dest,'GSE53697','LIBD',current)
            self.assertEqual((dest/'fractions.csv').read_bytes(),b'verified output')
            self.assertEqual(read_json(dest/'engine_reuse.json')['source_protocol'],digest(prior))
            with self.assertRaisesRegex(ValueError,'scientific protocol'):
                reuse_estimates(source,dest,'GSE53697','LIBD',dict(current,schema=3))
            changed=dict(current,code=dict(current['code'],**{'pipeline/clinical_estimators.R':'changed'}))
            with self.assertRaisesRegex(ValueError,'estimator/input code'):
                reuse_estimates(source,dest,'GSE53697','LIBD',changed)
            atomic_json(dest/'plan.json',dict(reference='anywhere',genes=['g1'],test_rows=[2]))
            with self.assertRaisesRegex(ValueError,'raw counts or estimator plan'):
                reuse_estimates(source,dest,'GSE53697','LIBD',current)
            atomic_json(dest/'plan.json',dict(reference='anywhere',genes=['g1'],test_rows=[1]))
            (dest/'raw.bin').write_bytes(b'changed raw')
            with self.assertRaisesRegex(ValueError,'raw counts or estimator plan'):
                reuse_estimates(source,dest,'GSE53697','LIBD',current)
            (dest/'raw.bin').write_bytes(b'raw')
            (old/'fractions.csv').write_bytes(b'corrupt')
            with self.assertRaises(ValueError):
                reuse_estimates(source,dest,'GSE53697','LIBD',current)


if __name__=='__main__':
    unittest.main()
