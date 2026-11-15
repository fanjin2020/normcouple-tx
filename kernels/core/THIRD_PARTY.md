# Third-party code and attribution

- MEAD 1.0.2: Dongyue Xie and contributors, GPL >=3. Unmodified source copied from https://github.com/DongyueXie/MEAD at commit `d66fa10491fb86de720a173e64b5cb26eeeefcd7`. The package DESCRIPTION omits imports; the installer supplies them separately. Native MEAD preprocessing is not claimed: this study freezes a reference-only, test-disjoint gene panel and uses the exported `MEAD_getX` and `MEAD_est` functions.
- vashr 0.99.1: Mengyin Lu and Matthew Stephens, GPL-3. Unmodified source archive from https://github.com/mengyin/vashr at commit `4eaae3edf86b3c0fadc3c57bcb65a71d2a9c3755`. The archive and all other vendored files are SHA-256 pinned in `bundle_manifest.json`.
- dtangle 2.0.9: official archive `https://raw.githubusercontent.com/gjhunt/dtangle/master/dtangle_2.0.9.tar.gz`; SHA-256 `ae665623fe98b463590f8abd3b5efbb768de28b195e7603998d517808e690129`. Copyright/license information is retained inside the archive.
- The unchanged Paper 3 stage-2 engine is the user's previous analysis source. SHA-256 `7ea5ea22bf76660f005d1c30ef40c5d0aaa8c1fffc1c0f24f3a67cd07b2f8011`.
- edgeR and DESeq2 are installed from Bioconductor and called through their official R functions, not reimplementations. See https://bioconductor.org/packages/edgeR and https://bioconductor.org/packages/DESeq2.

The new glue code is provided under GPL-3.0-or-later for compatibility with the vendored GPL components. Original copyrights remain with their respective authors. License text: https://www.gnu.org/licenses/gpl-3.0.html. No third-party research matrices are embedded in this code bundle.
