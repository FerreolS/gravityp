from __future__ import with_statement
from __future__ import absolute_import
from __future__ import print_function
import sys


try:
    import reflex
    from reflex_plot_widgets import *
    from gravity_plot_common import *
    import_success = True

except ImportError:
    import_success = False
    print("Error importing modules pyfits, wx, matplotlib, numpy")


def paragraph(text, width=None):
    """ wrap text string into paragraph
       text:  text to format, removes leading space and newlines
       width: if not None, wraps text, not recommended for tooltips as
              they are wrapped by wxWidgets by default
    """
    import textwrap
    if width is None:
        return textwrap.dedent(text).replace('\n', ' ').strip()
    else:
        return textwrap.fill(textwrap.dedent(text), width=width)


class DataPlotterManager(object):
    """
    This class must be added to the PipelineInteractiveApp with setPlotManager
    It must have following member functions which will be called by the app:
     - setInteractiveParameters(self)
     - readFitsData(self, fitsFiles):
     - addSubplots(self, figure):
     - plotProductsGraphics(self)
    Following members are optional:
     - setWindowHelp(self)
     - setWindowTitle(self)
     - setCurrentParameterHelper(self, helper)
    """

    # static members
    recipe_name = "gravi_viscal"

    def setInteractiveParameters(self):
        """
        This function specifies which are the parameters that should be presented
        in the window to be edited.  Note that the parameter has to also be in the
        in_sop port (otherwise it won't appear in the window). The descriptions are
        used to show a tooltip. They should match one to one with the parameter
        list.
        """
        return [
        ]

    def readFitsData(self, fitsFiles):
        """
        This function should be used to read and organize the raw fits files
        produced by the recipes.
        It receives as input a list of reflex.FitsFiles
        """
        # organize the files into a dictionary, here we assume we only have 
        # one file per category if there are more, one must use a
        # dictionary of lists
        self.frames = dict()
        calibrators = list()
        science = list()
        self.raw_tf = None
        for f in fitsFiles:
            if f.category == "SINGLE_CAL_VIS" or f.category == "DUAL_CAL_VIS":
                calibrators.append(f)
            elif f.category == "SINGLE_SCI_VIS" or f.category == "DUAL_SCI_VIS":
                science.append(f)
    
        if len(calibrators) != 0:
            self.raw_tf = PlotableRawCalibratorTF(calibrators)
    
        # we only have two states, we have data or we don't
        # define the plotting functions we want to use for each
        if self.raw_tf is not None: 
            self._add_subplots = self._add_subplots
            self._plot = self._data_plot
        else:
            self._add_subplots = self._add_nodata_subplots
            self._plot = self._nodata_plot

    def addSubplots(self, figure):
        """
        This function should be used to setup the subplots of the gui.  The the
        matplotlib documentation for a description of subplots.
        """
        self._add_subplots(figure)

    def plotProductsGraphics(self):
        """
        This function should be used to plot the data onto the subplots.
        """
        self._plot()

    def setWindowHelp(self):
        return 'Help for rrrecipe interactive window'

    def setWindowTitle(self):
        return 'gravity_viscal interactive window'

    def _add_nodata_subplots(self, figure):
        self.txt_plot = figure.add_subplot(111)

    def _add_subplots(self, figure):
        self.raw_tf_plots = list()
        self.raw_tf_plots.append(figure.add_subplot(321))
        self.raw_tf_plots.append(figure.add_subplot(322))
        self.raw_tf_plots.append(figure.add_subplot(323))
        self.raw_tf_plots.append(figure.add_subplot(324))
        self.raw_tf_plots.append(figure.add_subplot(325))
        self.raw_tf_plots.append(figure.add_subplot(326))

    def _data_plot(self):
        self.raw_tf.plot(self.raw_tf_plots, 
                         "Transfer function across the night",
                         "Transfer function across the night")
        

    def _nodata_plot(self):
        # could be moved to reflex library?
        self.txt_plot.set_axis_off()
        text_nodata = "Data not found. No calibrators found" 
        self.txt_plot.text(0.1, 0.6, text_nodata, color='#11557c',
                      fontsize=18, ha='left', va='center', alpha=1.0)
        self.txt_plot.tooltip = 'Data not found'

    def setCurrentParameterHelper(self, helper) :
        self.getCurrentParameterHelper = helper

#This is the 'main' function
if __name__ == '__main__':
    from reflex_interactive_app import PipelineInteractiveApp

    # Create interactive application
    interactive_app = PipelineInteractiveApp(enable_init_sop=True)

    # get inputs from the command line
    interactive_app.parse_args()

    #Check if import failed or not
    if not import_success:
        interactive_app.setEnableGUI(False)

    #Open the interactive window if enabled
    if interactive_app.isGUIEnabled():
        #Get the specific functions for this window
        dataPlotManager = DataPlotterManager()

        interactive_app.setPlotManager(dataPlotManager)
        interactive_app.showGUI()
    else:
        interactive_app.set_continue_mode()

    #Print outputs. This is parsed by the Reflex python actor to
    #get the results. Do not remove
    interactive_app.print_outputs()
    sys.exit()
