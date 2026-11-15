"""Real-R parity with the frozen bridge and independent exclusion controls."""
import os
from pathlib import Path
import sys
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'frozen'/'engine'),str(HERE/'frozen')]
import numpy as np
import pandas as pd
import revision_core as C
import revision_driver as D
from extension_driver import marker_properties
from selection import freeze
from reference_allocation import reference_stress
from experiment import exclusion_anchor,read_predictions,METHODS,INVARIANT
from common import read_json
from runtime import run_r


class BridgeParity(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('P3_TEST_RSCRIPT'),'Real R path required')
    def test_baseline_and_exclusion(self):
        with tempfile.TemporaryDirectory(prefix='p3_bridge_parity_') as tmp:
            root=Path(tmp)
            job=D.synthetic_jobs(root/'prepared')[0]
            freeze([job],1)
            job['budget_allocation']='reference'
            unit=C.plan_units([job],smoke=True)[0]
            base=C.draw_base(job,unit,0,C.SEEDS)
            data=C.intervene(job,base,dict(name='up025',effect=.25,pattern='up',panel='frozen'),C.SEEDS)
            stress,_=reference_stress(job,base,'legacy')
            A,_,anchor=marker_properties(job,base)
            axes=read_json(Path(job['reference_folder'])/'axes.json')
            raw=base['draw']['counts'][:,job['norm_genes']].T
            old=root/'old'; old.mkdir()
            for name,a in [('raw.bin',raw),('A.bin',A),('total.bin',data['original_total']),
                           ('add.bin',data['additions'].T),('stress.bin',stress.T)]:
                C.write_bin(old/name,a)
            entry=dict(key='r0',n=80,raw='raw.bin',A='A.bin',total='total.bin',marker_anchor=anchor,
                       outcomes=[dict(name='up025',additions='add.bin')],subsets=[dict(name='legacy',stress='stress.bin')])
            C.write_json(old/'jobs.json',dict(reference_folder=job['reference_folder'],entries=[entry]))
            exe,lib=os.environ['P3_TEST_RSCRIPT'],os.environ['P3_TEST_R_LIBRARY']
            os.environ['R_LIBS']=lib
            run_r([exe,'--vanilla',str(HERE/'frozen'/'bridge.R'),str(old),lib],lib,old/'R.log',120)
            oldrows=pd.read_csv(old/'outputs.csv',keep_default_na=False)
            new=root/'new'; new.mkdir()
            observed=raw.copy()
            observed[np.array(axes['marker_rows'])-1]+=stress.T
            observed[np.array(axes['test_rows'])-1]+=data['additions'].T
            for name,a in [('raw.bin',observed),('A.bin',A),('total.bin',data['total']+stress.sum(1))]:
                C.write_bin(new/name,a)
            steps=[0.,-.001,.001,-.0005,.0005,.01,.05,.10]
            entry.update(non_test_anchor=exclusion_anchor(job,base,A),group=base['group'].tolist(),steps=steps)
            C.write_json(new/'jobs.json',dict(axes=axes,entries=[entry]))
            run_r([exe,'--vanilla',str(HERE/'bridge.R'),str(new),lib],lib,new/'R.log',120)
            pp,_=read_predictions(new,{('r0',m,t) for m in METHODS for t in steps})
            for row in oldrows.itertuples():
                original=np.loadtxt(old/row.file,delimiter=',',ndmin=2)
                np.testing.assert_allclose(pp['r0',row.method,0.][0],original,atol=1e-12,rtol=1e-12)
            for method in INVARIANT:
                for t in steps:
                    np.testing.assert_array_equal(pp['r0',method,0.][0],pp['r0',method,t][0])


if __name__=='__main__':
    unittest.main()
