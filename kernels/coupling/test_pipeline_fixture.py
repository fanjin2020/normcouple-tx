"""Exercise orchestration with a labeled Python bridge double; never R validation."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'frozen'/'engine'),str(HERE/'frozen')]
import numpy as np
import pandas as pd
import revision_driver as D
import revision_core as C
import enhancement_experiment as E
from selection import freeze
from matched_markers import freeze_matched,variants
from enhancement_report import summarize
from common import save_json,file_manifest
from delivery import export_review


def bridge_double(command,library,log,timeout):
    root=Path(command[-2])
    plan=json.loads((root/'jobs.json').read_text())
    axes=plan['axes']
    mr=np.array(axes['marker_rows'])-1
    tr=np.array(axes['test_rows'])-1
    ng=axes['normalization_genes']
    rows=[]
    for entry in plan['entries']:
        n=entry['n']
        raw=np.fromfile(root/entry['raw'],dtype='<f8').reshape((ng,n),order='F')
        A=np.fromfile(root/entry['A'],dtype='<f8').reshape((len(mr),7),order='F')
        total=np.fromfile(root/entry['total'],dtype='<f8')
        group=np.array(entry['group'])
        marker=raw[mr].T
        for step in entry['steps']:
            library_=total-raw[tr].sum(0)+(raw[tr]*np.exp(step*group)[None,:]).sum(0)
            for method in E.DT_METHODS:
                if method=='tmm_exclude_both':
                    d=(total-raw[tr].sum(0))/entry['non_test_anchor']
                elif method=='marker_denominator_only':
                    d=marker.sum(1)/entry['marker_anchor']
                elif method=='marker_matched':
                    d=marker.sum(1)
                else:
                    d=library_
                ref=A/A.sum(0) if method=='marker_matched' else A
                q,reason=E.nnls_prediction(ref,marker,d)
                filename='%s_%s_%+.8f.csv'%(entry['key'],method,step)
                np.savetxt(root/filename,np.column_stack([d,q]),delimiter=',')
                rows.append(dict(key=entry['key'],method=method,step=step,file=filename,reason=reason,
                    warnings='PYTHON_BRIDGE_DOUBLE_NOT_DTANGLE',scale=1e6*(entry['marker_anchor'] if method=='marker_matched' else 1.)))
    pd.DataFrame(rows).to_csv(root/'outputs.csv',index=False)
    Path(log).write_text('Software fixture with Python NNLS adapter double; no R or TMM evidence.\n')


class PipelineFixture(unittest.TestCase):
    def test_batches_resume_reports_and_export(self):
        with tempfile.TemporaryDirectory(prefix='p3_review_fixture_') as tmp:
            root=Path(tmp)/'run'
            root.mkdir()
            jobs=D.synthetic_jobs(root/'prepared')
            freeze(jobs,1)
            freeze_matched(jobs,1,root/'matching')
            units=C.plan_units(jobs,smoke=True)[:1]
            settings=dict(C.SEEDS,repetitions=1,batch_size=1,random_subsets=1,
                subset_ids=variants(jobs[0]['stress_subsets']),diagnostic_subsets=['zero','legacy','random_001','matched_001'],
                derivative_step=.001,probe_steps=[.01],diagnostic_repetitions=1,software_fixture=True,phase='pilot')
            cfg=dict(signature='fixture_not_scientific_evidence',settings=settings,units=units)
            save_json(root/'protocol.json',cfg)
            E.JOBS=jobs
            E.OPTIONS=dict(settings,output=root,signature=cfg['signature'],r_info=dict(executable='NOT_R',library='NONE'),r_timeout=60)
            with patch.object(D,'run_r',bridge_double):
                first=E.run_batch((units[0],0))
                second=E.run_batch((units[0],0))
            self.assertFalse(first['resumed'])
            self.assertTrue(second['resumed'])
            result=summarize(root,cfg,'simulation')
            self.assertTrue(result['simulation_complete'])
            self.assertFalse(result['clinical_complete'])
            save_json(root/'final'/'STATUS.json',dict(status='SOFTWARE_FIXTURE_ONLY',**result))
            file_manifest(root/'final')
            self.assertEqual(export_review(root),0)
            self.assertTrue(list(root.parent.glob('p3r_*.zip')))
            summary=pd.read_csv(root/'final'/'method_summary.csv')
            self.assertEqual(set(summary.method),{'oracle',*E.METHODS})
            self.assertTrue((root/'final'/'matched_fdp.pdf').exists())
            marker_file=root/'final'/'marker_matching_balance.csv'
            marker_file.write_bytes(marker_file.read_bytes()+b'\n')
            with self.assertRaises(ValueError):
                export_review(root)


if __name__=='__main__':
    unittest.main()
