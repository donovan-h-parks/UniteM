###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import math
import logging

from unitem.utils import alphanumeric_sort

import svgwrite


class PlotCommonBases():
    """Create heatmap showing percentage of common bases between bins."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        # size of table labels
        self.font_size = 8

        # size of labels in gene cells
        self.cell_font_size = 6

        # width and height of each gene cell
        self.row_height = 2.5*self.font_size
        self.col_width = self.row_height
        self.cell_offset = 2  # border around each gene cell

        # expanded width of columns with comp and cont information
        self.qual_col_width = 4*self.font_size

    def _render_label_col(self, dwg,
                          labels,
                          start_x,
                          start_y,
                          group_id):
        """Render column of labels."""

        label_group = svgwrite.container.Group(id=group_id)
        dwg.add(label_group)

        for i, label in enumerate(labels):
            t = dwg.text(label,
                         x=[(start_x)],
                         y=[(start_y + i*self.row_height)],
                         font_size="%fpt" % self.font_size,
                         text_anchor='end',
                         direction='ltr',
                         fill='black')

            label_group.add(t)

    def _render_label_row(self, dwg,
                          labels,
                          label_start_x,
                          row_start_y,
                          group_id,
                          rotation=-45):
        """Render column of labels."""

        label_group = svgwrite.container.Group(id=group_id)
        dwg.add(label_group)

        for i, label in enumerate(labels):
            x = label_start_x + (i+0.5)*self.col_width
            t = dwg.text(label,
                         x=[(x)],
                         y=[(row_start_y)],
                         font_size="%fpt" % self.font_size,
                         text_anchor='start',
                         direction='ltr',
                         fill='black')

            t.rotate(rotation, (x, row_start_y))

            label_group.add(t)

    def _render_genome_quality_cols(self,
                                    dwg,
                                    bin_labels,
                                    quality,
                                    header_start_y,
                                    label_start_x,
                                    label_start_y):
        """Plot genome completeness and contamination columns."""

        # render header
        quality_group = svgwrite.container.Group(id='genome_quality')
        dwg.add(quality_group)

        for c, gene_name in enumerate(['Completeness (%)', 'Contamination (%)']):
            x = label_start_x + (c + 0.5)*self.qual_col_width
            t = dwg.text(gene_name,
                         x=[(x)],
                         y=[(header_start_y)],
                         font_size="%fpt" % self.font_size,
                         text_anchor='start',
                         direction='ltr',
                         fill='black')

            t.rotate(-45, (x, header_start_y))

            quality_group.add(t)

        # render rows
        for r, bin_label in enumerate(bin_labels):
            for c, q in enumerate(quality[bin_label]):
                t = dwg.text('%.1f' % q,
                             x=[(label_start_x + (c+0.5)*self.qual_col_width)],
                             y=[(label_start_y + r*self.row_height)],
                             font_size="%fpt" % self.font_size,
                             text_anchor='middle',
                             direction='ltr',
                             fill='black')
                quality_group.add(t)

    def _cell_properties(self, perc_common_bases):
        """Get desired color and size of heat map cell for a given percentage of common bases."""

        color = 'rgb(255,255,255)'
        size = self.col_width-2*self.cell_offset

        if perc_common_bases >= 100:
            color = 'rgb(165,15,21)'
        elif perc_common_bases >= 90:
            color = 'rgb(222,45,38)'
            size = math.sqrt(0.9)*size
        elif perc_common_bases >= 80:
            color = 'rgb(251,106,74)'
            size = math.sqrt(0.8)*size
        elif perc_common_bases >= 70:
            color = 'rgb(252,146,114)'
            size = math.sqrt(0.7)*size
        elif perc_common_bases >= 60:
            color = 'rgb(252,187,161)'
            size = math.sqrt(0.6)*size
        elif perc_common_bases >= 50:
            color = 'rgb(254,229,217)'
            size = math.sqrt(0.5)*size

        return color, size

    def _render_row(self,
                    dwg,
                    bin_labels,
                    bm_labels,
                    common_bases,
                    start_x,
                    start_y):
        """Render rows showing percentage of common bases."""

        table_group = svgwrite.container.Group(id='table')
        dwg.add(table_group)

        for r, bin_label in enumerate(bin_labels):
            for c, bm_label in enumerate(bm_labels):
                perc_cb = common_bases[bin_label].get(bm_label, 0)
                color, size = self._cell_properties(perc_cb)

                # plot heat map box
                x = start_x + c*self.col_width
                y = start_y + r*self.row_height

                if perc_cb > 0:
                    base_color, base_size = self._cell_properties(0)
                    rect = dwg.rect(insert=(x+0.5*(base_size - size),
                                            y+0.5*(base_size - size)),
                                    size=(size, size),
                                    fill=color)
                    rect.stroke(color='rgb(196,196,196)', width=0.1)
                    table_group.add(rect)

    def _render_legend(self, dwg):
        """Render legend."""

        legend_group = svgwrite.container.Group(id='legend')
        dwg.add(legend_group)

        x = self.fig_size_x
        y = 0
        _base_color, base_size = self._cell_properties(0)
        for perc_common_bases in [100, 90, 80, 70, 60, 50]:
            color, size = self._cell_properties(perc_common_bases)
            rect = dwg.rect(insert=(x+0.5*(base_size - size),
                                    y+0.5*(base_size - size)),
                            size=(size, size),
                            fill=color)
            rect.stroke(color='rgb(196,196,196)', width=0.1)
            legend_group.add(rect)

            legend_str = '%d' % perc_common_bases
            if perc_common_bases != 100:
                legend_str = '>' + legend_str

            t = dwg.text(legend_str,
                         x=[(x + base_size + self.cell_offset)],
                         y=[(y + 0.5*base_size + 0.5*self.font_size)],
                         font_size="%fpt" % self.font_size,
                         direction='ltr',
                         fill='rgb(0,0,0)')
            legend_group.add(t)

            y += size + self.cell_offset

    def plot(self, common_bases, quality, output_plot):
        """Create plot.

        Parameters
        ----------
        common_bases : d[unitem bid][binning method] -> percent common bases
          Percentage of common bases for each binning method.
        quality : d[unitem id] -> (completeness, contamination)
          Completeness and contamination of bins.
        output_plot : str
          Desired output file.
        """

        # sort bin labels
        bin_labels = alphanumeric_sort(common_bases)

        # sort binning method labels
        binning_methods = set()
        for bid in common_bases:
            for bm in common_bases[bid]:
                binning_methods.add(bm)
        bm_labels = alphanumeric_sort(binning_methods)

        # setup SVG drawing
        table_start_x = 0
        table_start_y = 0

        if not output_plot.endswith('.svg'):
            output_plot += '.svg'

        self.fig_size_x = table_start_x + len(bm_labels)*self.col_width
        self.fig_size_x += 2*self.qual_col_width + 0.5*self.col_width
        self.fig_size_y = table_start_y + len(bin_labels)*self.row_height

        dwg = svgwrite.Drawing(filename=output_plot,
                               size=(self.fig_size_x, self.fig_size_y),
                               profile='full',
                               style='font-family:Arial')
        dwg.set_desc(title='UniteM shared base pair plot.')

        # render legend
        self._render_legend(dwg)

        # render binning method labels
        label_start_x = table_start_x
        label_start_y = table_start_y + 0.5*self.row_height + 0.45*self.font_size
        self._render_label_col(dwg,
                               bin_labels,
                               label_start_x,
                               label_start_y,
                               'bin_labels')

        # render completeness and contamination
        header_start_y = table_start_y - 0.5*self.font_size
        label_start_x = table_start_x + 0.5*self.col_width
        label_start_y = table_start_y + 0.5*self.row_height + 0.45*self.font_size
        self._render_genome_quality_cols(dwg,
                                         bin_labels,
                                         quality,
                                         header_start_y,
                                         label_start_x,
                                         label_start_y)

        # move start of table to account for genome quality info
        table_start_x += 2*self.qual_col_width + 0.5*self.col_width

        # write header line
        label_start_x = table_start_x
        header_row_start_y = table_start_y - 0.5*self.font_size
        self._render_label_row(dwg,
                               bm_labels,
                               label_start_x,
                               header_row_start_y,
                               'binning_method_labels')

        # render gene complement row
        self._render_row(dwg,
                         bin_labels,
                         bm_labels,
                         common_bases,
                         table_start_x,
                         table_start_y)

        dwg.save()
