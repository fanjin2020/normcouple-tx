# Figure specification

Archived figure inputs are in `publication/data`; a new run supplies its own `final/data`. Tables and figures share CSVs through `catalog.json`. Plotting does not refit models, remove observations, smooth data or change uncertainty intervals.

## Typography and layout

- Double-column width: 6.85 inches (174 mm).
- Arial, with DejaVu Sans as a portable fallback.
- Axis labels: 9 pt; ticks and legends: 8 pt; lowercase panel labels: 11 pt bold.
- White background, minimal spines, light horizontal guides, aligned panel grids.
- Captions are kept outside the figure in `publication/figure_captions.json`.
- Manuscript text and tables are black on white. Color is reserved for figures.

## Visual encoding

| Series | Color | Encoding |
| --- | --- | --- |
| TMM | `#245792` | Circles |
| MRN | `#007F73` | Squares |
| Factor-only exclusion | `#705799` | Triangles |
| Complete exclusion | `#4E5862` | Squares |
| LIBD reference | Dark or blue | Solid lines, filled points or labeled columns |
| SEA-AD reference | Rose | Open points, different shapes/lines or labeled columns |

The equal-budget plot uses open symbols for zero-allowing allocation and filled symbols for positive-only allocation. The same method pair is joined within each subset. The ordered marker-dose heatmaps shade every tested cell and print values at amplitude 1.0. A small open circle denotes an absolute FDR point estimate above 5%, not a significance test. Sequential maps encode rates and diverging maps encode signed differences. Color ranges cover the displayed values. Marker-class scatter plots preserve every observation and the original deterministic horizontal jitter.

## Export and checks

```bash
python reproduce.py export --published --output ../paper3_figures --formats pdf png tiff
python -m unittest discover -s tests -p "test_*.py"
```

Each render checks panel labels, axis labels, annotations and legends for canvas overflow. PDF text uses embedded TrueType fonts. TIFF export is RGB, 600 dpi, LZW-compressed. The output `figure_sources.json` records dimensions, input files and captions for every figure.

The figure-source index in `pipeline/publication.py` is shared by the workbook and renderer. Fig. 6 displays the predefined subset, each random subset and the random mean on paired RMSE and FDR axes; the other three routes retain matched and unperturbed points. Fig. 7c shows signed TMM-to-MRN variance increments, not unsigned component shares. Fig. 8 uses a common hexbin density scale and reports the three categories that have discoveries in at least one fit. Fig. 9a,b use S20; Fig. 9c,d use S13. The renderer checks the shared full-family dose, family size and threshold counts before combining the two tables. Fig. S10 uses S18 and S19; Fig. S11 uses S21. The CBC command draws Fig. S12 from the paired assay records and held-out performance table (S22 and S23).

For a submission replay, use `--formats pdf png --workbook`. PNG is rendered at 300 dpi. PDF and SVG omit export dates, and SVG IDs use a fixed salt. Use the same software and fonts for byte-level comparisons; substituting a font changes text metrics.
