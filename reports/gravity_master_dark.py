from adari_core.plots.text import TextPlot
from adari_core.data_libs.master_dark_bias import MasterDarkBiasReport

from .gravity_util import GravitySetupInfo, extensions, data_reader

from .gravity_util import GravityReportMixin

import os


class GravityMasterDarkReport(GravityReportMixin, MasterDarkBiasReport):
    def __init__(self):
        super().__init__("gravity_master_dark")
        self.data_readers["master_im"] = data_reader

    def parse_sof(self):
        master_dark = None

        for filename, catg in self.inputs:
            if catg == "DARK" and master_dark is None:
                master_dark = filename

        file_lists = []
        if master_dark is not None:
            file_lists.append(
                {
                    "master_im": master_dark,
                }
            )
        return file_lists

    def generate_panels(self, **kwargs):
        vspace = 0.3
        panels = {}

        for ext1 in extensions:
            panels_1 = super().generate_raw_cuts_panels(
                master_im_ext=ext1,
                hist_clipping="val",
                hist_n_clipping=[-100, 100],
            )
            for i, (panel, panel_descr) in enumerate(panels_1.items()):
                master_im = self.hdus[i]["master_im"]
                fname = os.path.basename(str(master_im.filename()))
                fname_date = fname.removeprefix("GRAVI.").removesuffix(".fits")
                panel_descr["report_description"] = (
                    f"GRAVITY dark panel - "
                    f"{os.path.basename(panel_descr['master_im'])}, "
                    f"{panel_descr['master_im_ext']}"
                )
                panel_descr["report_name"] = f"GRAVITY_{fname_date}_{ext1}"

                t1 = TextPlot(columns=1, v_space=vspace)
                col1 = (
                    str(master_im["PRIMARY"].header.get("INSTRUME")),
                    "EXTNAME: " + ext1,
                    "PRO CATG: "
                    + str(master_im["PRIMARY"].header.get("HIERARCH ESO PRO CATG")),
                    "FILE NAME: " + fname,
                    "RAW1 NAME: "
                    + str(
                        master_im["PRIMARY"].header.get(
                            "HIERARCH ESO PRO REC1 RAW1 NAME"
                        )
                    ),
                )
                t1.add_data(col1)
                panel.assign_plot(t1, 0, 0, xext=2)

                t2 = TextPlot(columns=1, v_space=vspace, xext=1)
                col2 = GravitySetupInfo.dark(self.hdus[i]["master_im"])
                t2.add_data(col2)
                panel.assign_plot(t2, 2, 0, xext=1)

                if ext1 == "IMAGING_DATA_FT":
                    panel.pop(0, 2)

            panels.update(panels_1)

        return panels


rep = GravityMasterDarkReport()
