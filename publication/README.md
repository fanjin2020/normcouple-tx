# Publication sources

Use the root `reproduce.py` entry point. `export --work RUN` aggregates verified computations; `export --published` explicitly replays the archived article data.

The 21 supplementary tables and eleven additional figure or diagnostic datasets are stored once in `data/`. `data/catalog.json` maps plot inputs to CSVs; `table_schema.json` defines the tables. `published_manifest.json` records the producing framework, clinical probe provenance (S20), and the imaging inputs and matched support (S21).

`reproduce_figures.py`, `figure_panels.py`, `figure_data.py` and `figure_style.py` render Fig1-Fig9 and FigS1-FigS11. The renderer checks its actual data reads against the figure index in `pipeline/publication.py` and writes `figure_sources.json` with physical CSV names, views, captions and dimensions. The workbook's Guide uses this index even when no figures are drawn. `figure_captions.json` contains the matching full-data captions. A computed run without S20 labels Fig9 as finite-dose only; S11 requires imaging measurements (S21).

`mrn_tools.py` exposes the marker-restricted transformation, derivative, yield conversion and paired HC3 calculations tested by the software suite. These reusable functions are distinct from full refitting. Input panels, gene families, units and hashes are in `input_definitions/`.

Normalization, relative-count and unperturbed power experiments use 200 repetitions per unit. Common-support results use 50 repetitions and average three random and three matched subsets separately. Rates are fractions unless the column explicitly specifies percent (`_pct` in S20). Keep the outcome label `null` as a string when reading CSVs.
