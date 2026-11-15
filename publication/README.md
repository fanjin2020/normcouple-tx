# Publication sources

Use the root `reproduce.py` entry point. `export --work RUN` aggregates verified computations; `export --published` explicitly replays the archived article data.

`data/` contains one canonical copy of each numerical source: fifteen supplementary tables and eleven additional figure/text-check datasets. `data/catalog.json` defines named views used by the plotter, including two views of Table S7. `table_schema.json` contains table titles, notes and columns. `published_manifest.json` authenticates the archived data.

`reproduce_figures.py`, `figure_panels.py` and `figure_style.py` are plotting modules for Fig1-Fig9 and FigS1-FigS7. Every render writes `figure_sources.json` with physical CSV names, views, captions and dimensions. `figure_captions.json` contains the shared captions. Numerical rows are not duplicated in another JSON or workbook inside the source package.

`mrn_tools.py` exposes the marker-restricted transformation, derivative, yield conversion and paired HC3 calculations tested by the software suite. These reusable functions are distinct from full refitting. Input panels, gene families, units and hashes are in `input_definitions/`.

Normalization, relative-count and unperturbed power experiments use 200 repetitions per unit. Common-support results use 50 repetitions and separately average three random and three matched subsets. Rates are fractions. Read CSVs with `keep_default_na=False, na_values=['', 'NaN', 'nan']` so `null` remains a string.
