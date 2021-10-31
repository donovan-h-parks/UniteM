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

import logging

import svgwrite


class TreeCommonBases():
    """Create dendrogram showing common bases between bins."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        # size of labels
        self.font_size = 8
        self.scale_font_size = 6

        # width and height of each gene cell
        self.row_height = 1.3*self.font_size
        self.scale_interval_width = 15

        self.label_offset = 2

    def _render_scale(self, dwg):
        """Render scale."""

        scale_group = svgwrite.container.Group(id='scale')
        dwg.add(scale_group)

        scale_text_group = svgwrite.container.Group(id='scale_text')
        dwg.add(scale_text_group)

        # draw baseline
        scale_y_pos = self.fig_size_y + 0.5*self.row_height
        scale = dwg.path("M%f,%f" % (0, scale_y_pos + 0.25*self.row_height))
        scale.push("L%f,%f" % (0, scale_y_pos))
        scale.push("L%f,%f" % (self.fig_size_x, scale_y_pos))
        scale.push("L%f,%f" %
                   (self.fig_size_x, scale_y_pos + 0.25*self.row_height))
        scale.fill(color='none')
        scale.stroke(color='black', width=0.5)
        scale_group.add(scale)

        for s in [50, 60, 70, 80, 90, 100]:
            xpos = self.fig_size_x - ((100-s)/10)*self.scale_interval_width
            t = dwg.text(int(s),
                         x=[(xpos)],
                         y=[(scale_y_pos + 0.75*self.row_height + self.label_offset)],
                         font_size="%fpt" % self.scale_font_size,
                         text_anchor='middle',
                         fill='rgb(0,0,0)')
            scale_text_group.add(t)

            if s != 50 and s != 100:
                tick = dwg.line(start=(xpos, scale_y_pos),
                                end=(xpos, scale_y_pos + 0.25*self.row_height),
                                fill='black',
                                stroke_width=0.5)
                tick.stroke(color='black')
                scale_group.add(tick)

    def _render_labels(self, dwg, bin_id, common_bases):
        """Render labels."""

        label_group = svgwrite.container.Group(id='labels')
        dwg.add(label_group)

        t = dwg.text(bin_id,
                     x=[(self.fig_size_x + self.label_offset)],
                     y=[(self.fig_size_y + 0.5*self.font_size)],
                     font_size="%fpt" % self.font_size,
                     direction='ltr',
                     fill='rgb(0,0,0)')
        label_group.add(t)

        y = self.fig_size_y
        for bm, bid, cb in common_bases:
            y -= self.row_height
            t = dwg.text(bm,
                         x=[(self.fig_size_x + self.label_offset)],
                         y=[(y + 0.5*self.font_size)],
                         font_size="%fpt" % self.font_size,
                         direction='ltr',
                         fill='rgb(0,0,0)')
            label_group.add(t)

    def _render_tree(self, dwg, common_bases):
        """Render tree."""

        tree_group = svgwrite.container.Group(id='tree')
        dwg.add(tree_group)

        x_start = self.fig_size_x
        y_start = self.fig_size_y
        width_start = 0
        for r, (_bm, _bid, cb) in enumerate(common_bases):
            width = ((100-cb)/10)*self.scale_interval_width
            delta_width = width - width_start

            branch = dwg.path("M%f,%f" % (x_start, y_start))
            branch.push("L%f,%f" % (x_start-delta_width, y_start))
            branch.push("L%f,%f" % (x_start-delta_width,
                        self.fig_size_y - (r+1)*self.row_height))
            branch.push("L%f,%f" % (self.fig_size_x,
                        self.fig_size_y - (r+1)*self.row_height))
            branch.fill(color='none')
            branch.stroke(color='black', width=1)
            tree_group.add(branch)

            x_start = self.fig_size_x - width
            y_start = self.fig_size_y - r*self.row_height - 0.5*self.row_height
            width_start = width

        # render root branch
        root = dwg.line(start=(x_start, y_start),
                        end=(0, y_start),
                        fill='black',
                        stroke_width=1)
        root.stroke(color='black')
        tree_group.add(root)

    def plot(self, bin_id, common_bases, output_plot):
        """Create plot.

        Parameters
        ----------
        bin_id : str
          Identifier of bin compared to other bins.
        common_bases : d[binning method][bin id] -> percent common bases
          Percentage of common bases for each binning method.
        output_plot : str
          Desired output file.
        """

        # setup SVG drawing
        start_y = 0

        if not output_plot.endswith('.svg'):
            output_plot += '.svg'

        self.fig_size_x = 5*self.scale_interval_width  # 50, 60, 70, 80, 90, 100
        self.fig_size_y = start_y + len(common_bases)*self.row_height

        dwg = svgwrite.Drawing(filename=output_plot,
                               size=(self.fig_size_x, self.fig_size_y),
                               profile='full')
        dwg.set_desc(title='UniteM shared base pair tree.')

        self._render_scale(dwg)
        self._render_labels(dwg, bin_id, common_bases)
        self._render_tree(dwg, common_bases)

        dwg.save()
