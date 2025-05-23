from adari_core.data_libs.master_interfer_calib import MasterInterferCalibReport
from adari_core.utils.utils import fetch_kw_or_default
from adari_core.plots.text import TextPlot

from . import GravityReportMixin
from .gravity_util import GravitySetupInfo
import os


class GravityCalibratorReport(GravityReportMixin, MasterInterferCalibReport):
    def __init__(self):
        super().__init__("gravity_calibrator")

    def parse_sof(self):
        file_lists = []

        for filename, catg in sorted(self.inputs):
            if catg == "SINGLE_CAL_VIS" or catg == "DUAL_CAL_VIS":
                file_lists.append(
                   {
                     "oi_data": filename,
                   }
                )

        return file_lists

    def find_ext_index(self, hdu, ext):
        ind = []
        list_info = sorted(hdu.info(output=False), key=lambda item: item[2])
        for item in list_info:
            if item[1] == ext:
                ind.append(item[0])
        return ind

    def generate_panels(self, **kwargs):

        extensions = ["OI_VIS2", "OI_T3", "OI_VIS", "OI_FLUX", "TF2"]
        # Wavelength
        wave_iext = self.find_ext_index(self.hdus[0]["oi_data"], "OI_WAVELENGTH")

        vspace = 0.4

        # meta text
        meta0 = []
        meta1 = []
        
        for ext in extensions:
            cols0 = []
            cols1 = []
            
            for hdu in self.hdus:
                mdata_hdul = hdu["oi_data"]
                mdata_ext = 0
                oi_data_iext = self.find_ext_index(hdu["oi_data"], ext)
                col0 = (
                    "INSTRUME: " + str(mdata_hdul[mdata_ext].header.get("INSTRUME")),
                    "EXTNAME: " + ext,
                    "PRO.CATG: " + str(fetch_kw_or_default(mdata_hdul[mdata_ext], "HIERARCH ESO PRO CATG", default="N/A")),
                    "FILE NAME: " + os.path.basename(mdata_hdul.filename()), 
                    "RAW1 NAME: " + str(fetch_kw_or_default(mdata_hdul[mdata_ext], "HIERARCH ESO PRO REC1 RAW1 NAME", default="N/A")),
                )
                for i in oi_data_iext:
                    col0insname = col0 + ("INSNAME: " + str(fetch_kw_or_default(mdata_hdul[i], "INSNAME", default="N/A")),)
                    cols0.append(col0insname)
                    cols1.append(GravitySetupInfo.calib(mdata_hdul))
            meta0.append(cols0)
            meta1.append(cols1)

        # 1. Visibility
        oi_data_iext = self.find_ext_index(self.hdus[0]["oi_data"], extensions[0])

        panel_vis2 = super().generate_panels(
           oi_ext_list=oi_data_iext,
           wave_ext_list=wave_iext,
           y_label="Visibility",
           y_column="VIS2DATA",
           y_err="VIS2ERR",
           y_lim="clip0",
           legend_from=True,
           stretch=0.8,
           **kwargs
           )
        for i, (panel, panel_descr) in enumerate(panel_vis2.items()):
            t0 = TextPlot(columns=1, v_space=vspace)
            t0.add_data(meta0[0][i])
            panel.assign_plot(t0, 0, 0, xext=1)

            t1 = TextPlot(columns=1, v_space=vspace, xext=1)
            t1.add_data(meta1[0][i][:6])
            panel.assign_plot(t1, 3, 0, xext=1)
        
            t2 = TextPlot(columns=1, v_space=vspace, xext=1)
            t2.add_data(meta1[0][i][6:])
            panel.assign_plot(t2, 4, 0, xext=1)

        # 2. Flux (if OI_FLUX)
        oi_data_iext = self.find_ext_index(self.hdus[0]["oi_data"], extensions[3])
        
        # check if OI_FLUX exists
        panel_oiflux = {}
        if len(oi_data_iext) != 0:
            panel_oiflux = super().generate_panels(
                oi_ext_list=oi_data_iext,
                wave_ext_list=wave_iext,
                n_panels=4,
                y_label="Flux",
                y_column="FLUX",
                y_err="FLUXERR",
                y_lim="clip0",
                legend_from=True,
                stretch=0.8,
                **kwargs
                )
            for i, (panel, panel_descr) in enumerate(panel_oiflux.items()):
                t0 = TextPlot(columns=1, v_space=vspace)
                t0.add_data(meta0[3][i])
                panel.assign_plot(t0, 0, 0, xext=1)

                t1 = TextPlot(columns=1, v_space=vspace, xext=1)
                t1.add_data(meta1[3][i][:6])
                panel.assign_plot(t1, 3, 0, xext=1)

                t2 = TextPlot(columns=1, v_space=vspace, xext=1)
                t2.add_data(meta1[3][i][6:])
                panel.assign_plot(t2, 4, 0, xext=1)

    
        # 3. Transfer function (if TF2)
        oi_data_iext = self.find_ext_index(self.hdus[0]["oi_data"], extensions[4])
        
        # check if TF2 exists
        panel_tf2 = {}
        if len(oi_data_iext) != 0:
            panel_tf2 = super().generate_panels(
                oi_ext_list=oi_data_iext,
                wave_ext_list=wave_iext,
                y_label="TF2",
                y_column="TF2",
                y_err="TF2ERR",
                y_lim="clip0",
                legend_from=True,
                stretch=0.8,
                **kwargs
                )
            for i, (panel, panel_descr) in enumerate(panel_tf2.items()):
                t0 = TextPlot(columns=1, v_space=vspace)
                t0.add_data(meta0[4][i])
                panel.assign_plot(t0, 0, 0, xext=1)

                t1 = TextPlot(columns=1, v_space=vspace, xext=1)
                t1.add_data(meta1[4][i][:6])
                panel.assign_plot(t1, 3, 0, xext=1)

                t2 = TextPlot(columns=1, v_space=vspace, xext=1)
                t2.add_data(meta1[4][i][6:])
                panel.assign_plot(t2, 2, 0, xext=1)


        panels = {**panel_vis2}
        if panel_oiflux:
            panels.update(panel_oiflux)
        if panel_tf2:
            panels.update(panel_tf2)

        return panels

rep = GravityCalibratorReport() 
