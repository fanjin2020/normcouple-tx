# Third-party attribution

The publication and orchestration code uses NumPy, pandas, SciPy, h5py and Matplotlib. Optional Excel export uses openpyxl. The R analyses call official dtangle, edgeR, DESeq2, and MEAD functions. Their package licenses remain applicable.

- MEAD 1.0.2, Dongyue Xie and contributors, GPL >= 3; commit `d66fa10491fb86de720a173e64b5cb26eeeefcd7`, https://github.com/DongyueXie/MEAD.
- vashr 0.99.1, Mengyin Lu and Matthew Stephens, GPL-3; commit `4eaae3edf86b3c0fadc3c57bcb65a71d2a9c3755`, https://github.com/mengyin/vashr.
- dtangle 2.0.9, official source archive retained with its package license: https://raw.githubusercontent.com/gjhunt/dtangle/master/dtangle_2.0.9.tar.gz. SHA256: `ae665623fe98b463590f8abd3b5efbb768de28b195e7603998d517808e690129`.
- The study's frozen kernels are retained byte-for-byte, including embedded source and input definitions. Their third-party notices remain in `kernels/`; `kernels/layout.json` records original internal layouts.

The glue-code notice supplied with the scientific sources specifies GPL-3.0-or-later. This package preserves that license instead of applying a conflicting permissive license to the bundled sources. Dataset and annotation permissions remain separate from the code license.
