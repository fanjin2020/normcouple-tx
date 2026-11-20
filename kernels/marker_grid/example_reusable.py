"""Small software example; replace these synthetic inputs for research use."""
import numpy as np
from mrn_tools import mrn_inputs, yield_convert, actual_span_hc3

rng=np.random.default_rng(14)
raw=rng.poisson(25,size=(80,300))
reference=rng.uniform(.00001,.0002,size=(300,7))
yields=np.linspace(2000,8000,7)
inputs=mrn_inputs(raw,reference,yields,anchor=.03)
print('MRN target and reference shapes:',inputs['target_log2'].shape,inputs['reference_log2'].shape)
# In real use, supply the official deconvolution output for these matched inputs.
rna=rng.dirichlet(np.ones(7),80)
props=yield_convert(rna,inputs['marker_yield'])
x=np.repeat([0.,1.],40)
y=.2*x[:,None]+rng.normal(size=(80,5))
fit=actual_span_hc3(y,x,props)
print('SOFTWARE EXAMPLE coefficient:',fit['beta'])
print('Available:',fit['available'],'nuisance rank:',fit['nuisance_rank'])
