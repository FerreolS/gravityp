from adari_core.data_libs.master_rawdisp import MasterRawdispReport
import os
from .gravity_util import GravitySetupInfo, data_reader
from adari_core.plots.cut import CutPlot
from adari_core.plots.images import ImagePlot
from adari_core.utils.utils import fetch_kw_or_default
from . import GravityReportMixin


class GravityRawdispReport(GravityReportMixin, MasterRawdispReport):
    def __init__(self):
        super().__init__("gravity_rawdisp")
        self.extensions = []
        self.tasks = {
            "DARK_RAW": "dark",
            #
            "FLAT_RAW": "pixel_to_visibility",
            "WAVESC_RAW": "pixel_to_visibility",
            "WAVE_RAW": "pixel_to_visibility",
            "P2VM_RAW": "pixel_to_visibility",
            #
            "SINGLE_CAL_RAW": "calibrator_visibilities",
            "DUAL_CAL_RAW": "calibrator_visibilities",
            "SINGLE_SKY_RAW": "calibrator_visibilities",
            "DUAL_SKY_RAW": "calibrator_visibilities",
            #
            "WAVELAMP_RAW": "arc",
            #
            "PIEZOTF_RAW": "actuators_response",
        }
        self.data_readers["filename"] = data_reader
        self.setup_info = GravitySetupInfo

    def parse_sof(self):
        # we building multiple report sets, so we append multiple reports to file_lists
        # get a list of tags
        ext = "IMAGING_DATA_SC"
        tags = list(self.tasks.keys())
        added = {}
        file_lists = []
        for filename, catg in self.inputs:
            if catg in tags:
                if filename is not None and catg not in added:
                    file_lists.append({"filename": filename})
                    added[catg] = self.tasks[catg]
                    self.sof_tag.append(catg)
                    self.extensions.append(ext)        
        return file_lists


    def generate_panels(self, **kwargs):
        panels = {}
        new_panels = super().generate_panels(ext=self.extensions, **kwargs)
        for i, (panel, panel_descr) in enumerate(new_panels.items()):
            # Alter the cut pos, or remove CutPlot(s) completely,
            # depending on task name
            try:
                task_name = panel_descr["task_name"]
                if task_name != "dark":
                    panel.pop(8,3) # delete X-direction cut
            except KeyError as e:
                raise RuntimeError(
                    "A report has been created by "
                    "MasterRawdispReport that did "
                    "not come back with a task name "
                    "attached!"
                )
            panel_descr["report_name"] = "gravity_rawdisp_{}_{}_{}_{}".format(
                task_name,
                self.extensions[i],
                self.sof_tag[i].lower(),
                os.path.basename(panel_descr["filename"]),
            )
            panel_descr["report_description"] = (
                f"GRAVITY rawdisp panel - "
                f"{panel_descr['task_name']}, "
                f"{panel_descr['tag']}, "
                f"{os.path.basename(panel_descr['filename'])}, "
                f"{panel_descr['ext']}"
            )
        panels = {**panels, **new_panels}

        return panels


rep = GravityRawdispReport()

