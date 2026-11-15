from pathlib import Path
import sys
import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'publication'))
import figure_style as style


class FigureStyleTests(unittest.TestCase):
    def tearDown(self):
        plt.close('all')

    def test_method_colors_are_consistent(self):
        self.assertEqual(style.METHOD_COLORS['tmm'],style.METHOD_COLORS['TMM'])
        self.assertNotEqual(style.METHOD_COLORS['tmm'],style.METHOD_COLORS['marker_matched'])
        self.assertEqual(len(set(style.REFERENCE_COLORS.values())),2)

    def test_dark_text_contrast(self):
        rgb=to_rgb(style.INK)
        linear=[v/12.92 if v<=.04045 else ((v+.055)/1.055)**2.4 for v in rgb]
        luminance=sum(a*b for a,b in zip(linear,[.2126,.7152,.0722]))
        self.assertGreater(1.05/(luminance+.05),7)

    def test_wide_panel_label_is_not_clipped(self):
        fig=plt.figure(figsize=(6.85,4))
        ax=fig.add_axes([.12,.2,.86,.6])
        style.panel(ax,'a')
        style.check_text_bounds(fig)
        box=ax.texts[0].get_window_extent(fig.canvas.get_renderer())
        self.assertGreaterEqual(box.x0,0)

    def test_outside_annotation_is_rejected(self):
        fig,ax=plt.subplots()
        ax.text(-2,1,'clipped',transform=ax.transAxes)
        with self.assertRaisesRegex(ValueError,'Text outside figure canvas'):
            style.check_text_bounds(fig)


if __name__=='__main__':
    unittest.main()
