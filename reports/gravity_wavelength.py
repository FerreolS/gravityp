from adari_core.report import AdariReportBase
from adari_core.plots.panel import Panel
from adari_core.plots.points import LinePlot
from adari_core.plots.text import TextPlot
from .gravity_util import GravitySetupInfo

import os

from .gravity_util import GravityReportMixin


class GravityArcReport(GravityReportMixin, AdariReportBase):
    def __init__(self):
        super().__init__("gravity_wavelength")

    def parse_sof(self):
        wavelamp_raw = None
        wavelamp = None

        for filename, catg in self.inputs:
            if catg == "WAVELAMP" and wavelamp is None:
                wavelamp = filename
            elif catg == "WAVELAMP_RAW" and wavelamp_raw is None:
                wavelamp_raw = filename

        # Build and return the file name list
        file_lists = []
        if wavelamp_raw is not None:
            if wavelamp is not None:
                file_lists.append(
                    {
                        "wavelamp": wavelamp,
                        "wavelamp_raw": wavelamp_raw,
                    }
                )

        return file_lists

    def plot_spectra(self, n, label):
        dx_arr = []
        for i in self.argon_pos:
            dx = [i, i]
            dx_arr.append(dx)

        dy_arr = []
        for j in self.argon_pos:
            dy = [
                0.5 * max(list(self.spectrum_data[0][n])),
                0.75 * max(list(self.spectrum_data[0][n])),
            ]
            dy_arr.append(dy)
        plot = LinePlot(legend=False, title=label)
        plot.add_data(
            d=[self.wave1, list(self.spectrum_data[0][n])], label=label, color="blue"
        )
        for k in range(len(dx_arr)):
            plot.add_data(
                d=[dx_arr[k], dy_arr[k]], label="Argon lines" + str(k), color="red"
            )

        return plot

    def first_panel(self):
        p = Panel(6, 5, height_ratios=[1, 4, 4, 4, 4])

        dx_arr = []
        for i in self.argon_pos:
            dx = [i, i]
            dx_arr.append(dx)

        dy_arr = []
        for j in self.argon_pos:
            dy = [
                0.5 * max(list(self.spectrum_data[0][0])),
                0.75 * max(list(self.spectrum_data[0][0])),
            ]
            dy_arr.append(dy)

        plot = LinePlot(legend=False, y_label="ADU", title="DATA1")
        plot.add_data(
            d=[self.wave1, list(self.spectrum_data[0][0])], label="DATA1", color="blue"
        )
        for k in range(len(dx_arr)):
            plot.add_data(
                d=[dx_arr[k], dy_arr[k]], label="Argon lines" + str(k), color="red"
            )

        p.assign_plot(plot, 0, 1)

        for i in range(1, 6):
            plot = self.plot_spectra(2 * i, "DATA" + str(i + 1))
            p.assign_plot(plot, i, 1)

        for i in range(0, 6):
            plot = self.plot_spectra((2 * i) + 6, "DATA" + str(i + 7))
            p.assign_plot(plot, i, 2)

        for i in range(0, 6):
            plot = self.plot_spectra((2 * i) + 12, "DATA" + str(i + 13))
            p.assign_plot(plot, i, 3)

        for i in range(0, 6):
            dy_arr = []
            for j in self.argon_pos:
                dy = [
                    0.5 * max(list(self.spectrum_data[0][(2 * i) + 18])),
                    0.75 * max(list(self.spectrum_data[0][(2 * i) + 18])),
                ]
                dy_arr.append(dy)

            plot = LinePlot(
                legend=False, x_label="Lambda [um]", title="DATA" + str(i + 19)
            )
            plot.add_data(
                d=[self.wave1, list(self.spectrum_data[0][(2 * i) + 18])],
                label="DATA" + str(i + 19),
                color="blue",
            )
            for k in range(len(dx_arr)):
                plot.add_data(
                    d=[dx_arr[k], dy_arr[k]], label="Argon lines" + str(k), color="red"
                )
            p.assign_plot(plot, i, 4)

        return p

    def second_panel(self):
        p = Panel(6, 5, height_ratios=[1, 4, 4, 4, 4])

        dx_arr = []
        for i in self.argon_pos:
            dx = [i, i]
            dx_arr.append(dx)

        dy_arr = []
        for j in self.argon_pos:
            dy = [
                0.5 * max(list(self.spectrum_data[0][24])),
                0.75 * max(list(self.spectrum_data[0][24])),
            ]
            dy_arr.append(dy)

        plot = LinePlot(
            legend=False,
            # x_label="Lambda [nm]",
            y_label="ADU",
            title="DATA25",
        )
        plot.add_data(
            d=[self.wave2, list(self.spectrum_data[0][24])],
            label="DATA25",
            color="blue",
        )
        for k in range(len(dx_arr)):
            plot.add_data(
                d=[dx_arr[k], dy_arr[k]], label="Argon lines" + str(k), color="red"
            )
        p.assign_plot(plot, 0, 1)

        for i in range(1, 6):
            plot = self.plot_spectra((2 * i) + 24, "DATA" + str(i + 25))
            p.assign_plot(plot, i, 1)

        for i in range(0, 6):
            plot = self.plot_spectra((2 * i) + 30, "DATA" + str(i + 31))
            p.assign_plot(plot, i, 2)

        for i in range(0, 6):
            plot = self.plot_spectra((2 * i) + 36, "DATA" + str(i + 37))
            p.assign_plot(plot, i, 3)

        for i in range(0, 6):
            dy_arr = []
            for j in self.argon_pos:
                dy = [
                    0.5 * max(list(self.spectrum_data[0][(2 * i) + 42])),
                    0.75 * max(list(self.spectrum_data[0][(2 * i) + 42])),
                ]
                dy_arr.append(dy)

            plot = LinePlot(
                legend=False,
                x_label="Lambda [um]",
                # y_label="ADU",
                title="DATA" + str(i + 43),
            )
            plot.add_data(
                d=[self.wave2, list(self.spectrum_data[0][(2 * i) + 42])],
                label="DATA" + str(i + 43),
                color="blue",
            )
            for k in range(len(dx_arr)):
                plot.add_data(
                    d=[dx_arr[k], dy_arr[k]], label="Argon lines" + str(k), color="red"
                )
            p.assign_plot(plot, i, 4)

        return p

    def generate_panels(self, **kwargs):
        panels = {}
        wavelamp = self.hdus[0]["wavelamp"]
        ext_data = "SPECTRUM_DATA_SC"
        self.spectrum_data = wavelamp[ext_data].data
        self.wave1 = wavelamp[4].data["EFF_WAVE"] / 1e-6
        self.wave2 = wavelamp[5].data["EFF_WAVE"] / 1e-6
        self.argon_pos = wavelamp["POS_ARGON"].data["WAVE_TH"] / 1e-6
        fname = wavelamp.filename()
        hdr = wavelamp[0].header
        hdr_oi_1 = wavelamp[4].header
        hdr_oi_2 = wavelamp[5].header
        wavelamp_procatg = hdr.get("HIERARCH ESO PRO CATG")

        px = 0
        vspace = 0.2

        t1a = TextPlot(columns=1, v_space=vspace)
        col1a = (
            str(hdr.get("INSTRUME")),
            "EXTNAME: " + ext_data,
            "PRO CATG: " + str(hdr.get("HIERARCH ESO PRO CATG")),
            "FILE NAME: " + os.path.basename(fname),
            "RAW1.NAME: " + str(hdr.get("HIERARCH ESO PRO REC1 RAW1 NAME")),
            "INSNAME: " + str(hdr_oi_1.get("INSNAME")),
        )
        t1a.add_data(col1a)

        t1b = TextPlot(columns=1, v_space=vspace)
        col1b = (
            str(hdr.get("INSTRUME")),
            "EXTNAME: " + ext_data,
            "PRO CATG: " + str(hdr.get("HIERARCH ESO PRO CATG")),
            "FILE NAME: " + os.path.basename(fname),
            "RAW1.NAME: " + str(hdr.get("HIERARCH ESO PRO REC1 RAW1 NAME")),
            "INSNAME: " + str(hdr_oi_2.get("INSNAME")),
        )
        t1b.add_data(col1b)

        px = px + 2
        t2 = TextPlot(columns=1, v_space=vspace, xext=1)
        self.metadata = GravitySetupInfo.arc(wavelamp)
        col2 = self.metadata[:3]
        t2.add_data(col2)

        px = px + 2
        t3 = TextPlot(columns=1, v_space=vspace, xext=1)
        col3 = self.metadata[3:]
        t3.add_data(col3)

        pola_mode = hdr.get("HIERARCH ESO FT POLA MODE")

        input_files = [wavelamp.filename()]
        if pola_mode == "SPLIT":
            p1 = self.first_panel()
            p1.assign_plot(t1a, 0, 0, xext=2)
            p1.assign_plot(t2, 2, 0)
            p1.assign_plot(t3, 4, 0)

            p2 = self.second_panel()
            p2.assign_plot(t1b, 0, 0, xext=2)
            p2.assign_plot(t2, 2, 0)
            p2.assign_plot(t3, 4, 0)

            addme1 = {
                "panel": "first panel",
                "report_name": f"gravity_{wavelamp_procatg}_1",
                "report_description": "GRAVITY arc wavelength panel 1",
                "report_tags": [],
                "input_files": input_files                
            }

            addme2 = {
                "panel": "second panel",
                "report_name": f"gravity_{wavelamp_procatg}_2",
                "report_description": "GRAVITY arc wavelength panel 2",
                "report_tags": [],
                "input_files": input_files
            }

            panels[p1] = addme1

            panels[p2] = addme2

        else:
            p1 = self.first_panel()
            p1.assign_plot(t1a, 0, 0, xext=2)
            p1.assign_plot(t2, 2, 0)
            p1.assign_plot(t3, 4, 0)

            addme = {
                "report_name": f"gravity_{wavelamp_procatg}",
                "report_description": "GRAVITY arc wavelength panel",
                "report_tags": [],
                "input_files": input_files
            }

            panels[p1] = addme

        return panels


rep = GravityArcReport()
