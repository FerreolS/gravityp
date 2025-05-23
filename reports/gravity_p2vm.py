from adari_core.report import AdariReportBase
from adari_core.plots.panel import Panel
from adari_core.plots.points import ScatterPlot, LinePlot
from adari_core.plots.histogram import HistogramPlot
from adari_core.plots.cut import CutPlot
from adari_core.plots.text import TextPlot
from adari_core.utils.utils import fetch_kw_or_default
from adari_core.plots.images import ImagePlot
from .gravity_util import GravitySetupInfo
from astropy.io import fits


import numpy as np
import os

from . import GravityReportMixin


class GravityP2vmReport(GravityReportMixin, AdariReportBase):
        
    def __init__(self):
        super().__init__("gravity_p2vm")

    def parse_sof(self):
        files_list = []
        
        result = {
            "flat": None,
            "p2vm": None
        }
        flat_raw_counter = 1
        p2vm_raw_counter = 1
        
        for filename, catg in sorted(self.inputs):
            if catg == "FLAT":
                result["flat"] = filename
            elif catg == "P2VM":
                result["p2vm"] = filename
            elif catg == "FLAT_RAW":
                result[f"flat_raw_{flat_raw_counter}"] = filename
                flat_raw_counter += 1
            elif catg == "P2VM_RAW":
                result[f"p2vm_raw_{p2vm_raw_counter}"] = filename
                p2vm_raw_counter += 1
                
        file_lists = [result]
        return file_lists

    def panel_images_cuts_sc(self, n, prod_name, include_product):
        """Function to generate image plots, cut plots, and histogram plots for SC, 
        used to generate panels 1 to 4.

        Inputs: 
        n - number of rows
        raw_name - key name of raw file in files list dict
        prod_name - key name of product file in files list dict
        include_product - 1 if product file is included in panel, 0 if not.
        
        Outputs:
        Panel objects p1 (image plots) and p3 (cut plots and histograms)
        """
        ext_data_sc = "IMAGING_DATA_SC"
        ext_detector_sc = "IMAGING_DETECTOR_SC"
        ext_oi_sc = "OI_WAVELENGTH"

        p2vm = self.hdus[0]["p2vm"]
        self.metadata = GravitySetupInfo.p2vm(p2vm)
        hdr = self.hdus[0][prod_name]["PRIMARY"].header
        hdr_sc = self.hdus[0][f"{prod_name}_raw_1"][ext_data_sc].header
        self.fname = self.hdus[0][prod_name].filename()
        self.procatg = hdr.get("HIERARCH ESO PRO CATG")


        #Initiate panels
        p1 = Panel(n, 2, height_ratios=[1, 4])

        p3 = Panel(4, 7, height_ratios=[1, 4, 4, 4, 4, 4, 4])

        #Text plots
        vspace = 0.2
        detector_list = list(self.hdus[0][f"{prod_name}_raw_1"][ext_detector_sc].data)
        dname_list = []
        dpos_list = []
        for i in range(len(detector_list)):
            dname = detector_list[i][4]
            dpos = detector_list[i][5][1]
            dname_list.append(dname)
            dpos_list.append(dpos)

        t1a = TextPlot(columns=1, v_space=vspace)
        col1a = (
            str(hdr.get("INSTRUME")),
            "EXTNAME: " + ext_data_sc, 
            "PRO CATG: "
            + str(hdr.get("HIERARCH ESO PRO CATG")),
            "FILE NAME: " + os.path.basename(self.fname),
            "RAW1.NAME: " + str(hdr.get("HIERARCH ESO PRO REC1 RAW1 NAME")),
            "NAXIS3: " + str(hdr_sc.get("NAXIS3")),
        )
        t1a.add_data(col1a)
        p1.assign_plot(t1a, 0, 0, xext=2)
        p3.assign_plot(t1a, 0, 0, xext=2)

        t2 = TextPlot(columns=1, v_space=vspace, xext=1)
        col2 = self.metadata[:4]
        t2.add_data(col2)
        p1.assign_plot(t2, 2, 0, xext=1)
        p3.assign_plot(t2, 2, 0, xext=1)

        t3 = TextPlot(columns=1, v_space=vspace, xext=1)
        col3 = self.metadata[4:]
        t3.add_data(col3)
        p1.assign_plot(t3, 3, 0)
        p3.assign_plot(t3, 3, 0)

        if prod_name == "p2vm":
            nfiles = 6
        else:
            nfiles = 4

        #Raw image plots
        for i in range(nfiles):
            im_data = self.hdus[0][f"{prod_name}_raw_{i+1}"][ext_data_sc].data
            a = np.shape(im_data.mean(axis=0))[1] / np.shape(im_data.mean(axis=0))[0]
            #Generate image plot
            im_plot= ImagePlot(
                im_data.mean(axis=0), 
                title=f"{prod_name.upper()}_RAW_{str(i+1)}",
                v_clip="percentile",
                v_clip_kwargs={"percentile": 95.0},
                #show_colorbar=False,
                aspect=2*a
            )
            p1.assign_plot(im_plot, i, 1)

            #Generate cut plot
            cutpos = im_plot.get_data_coord(im_plot.data.shape[1] // 2, "x")
            im_cut = CutPlot(
                "x",
                title=f"columns (raw #{i+1}@ X {cutpos})",
                legend = False, 
                y_label=fetch_kw_or_default(
                    self.hdus[0][prod_name]["PRIMARY"], "BUNIT", default="ADU"
                )
            )
            im_cut.add_data(im_plot, cutpos, color="black", label="master")
            p3.assign_plot(im_cut, 0, i+1, xext=3)\

            #Generate histogram
            raw_hist = HistogramPlot(
                title=f"raw #{i+1}",
                x_label=fetch_kw_or_default(self.hdus[0][prod_name]["PRIMARY"], "BUNIT", "counts"),
                v_clip = "val",
                v_clip_kwargs = {"low": -1000, "high": 30000}
            )
            raw_hist.add_data(im_plot, color="black")
            p3.assign_plot(raw_hist, 3, i+1)

        #Generate optional product plots
        if include_product == 1:
    
            nx = np.shape(self.hdus[0][prod_name][ext_data_sc].data[0,:,:])[0]
            ny = np.shape(self.hdus[0][prod_name][ext_data_sc].data[0,:,:])[1]
    
            im_plot2 = ImagePlot(
                self.hdus[0][prod_name][ext_data_sc].data[0,:,:], 
                title=f"{self.procatg}",
                v_clip="percentile",
                v_clip_kwargs={"percentile": 95.0},
                y_ticks=dpos_list,
                aspect = 2*(ny / nx)
                )
    
            p1.assign_plot(im_plot2, n-1, 1)

            cutpos = im_plot2.get_data_coord(im_plot2.data.shape[1] // 2, "x")
                
            im_cut = CutPlot(
                "x",
                title="columns (master @ X {})".format(cutpos),
                legend = False,
                y_label=fetch_kw_or_default(
                    self.hdus[0][prod_name]["PRIMARY"], "BUNIT", default="ADU"
                )
            )
            im_cut.add_data(im_plot2, cutpos, color="red", label="master")
            p3.assign_plot(im_cut, 0, n, xext=3)
            
            master_hist = HistogramPlot(
                title="flat",
                x_label=fetch_kw_or_default(self.hdus[0][prod_name]["PRIMARY"], "BUNIT", "counts"),
                color="red",
                v_clip = "val",
                v_clip_kwargs = {"low": -1000, "high": 30000}
            )
            master_hist.add_data(im_plot2, color="red")
            p3.assign_plot(master_hist, 3, n)


        return p1, p3


    def panel_images_cuts_ft(self, n, prod_name, include_product):
        """Function to generate image plots for FT, 
        used to generate panels 5 and 6.

        Inputs: 
        n - number of rows
        raw_name - key name of raw file in files list dict
        prod_name - key name of product file in files list dict
        include_product - 1 if product file is included in panel, 0 if not.
        
        Outputs:
        Panel object p2 (image plots)
        """

        ext_data_ft = "IMAGING_DATA_FT"
        ext_detector_ft = "IMAGING_DETECTOR_FT"        
        #self.flat = self.hdus[0]["flat"]
        p2vm = self.hdus[0]["p2vm"]
        self.metadata = GravitySetupInfo.p2vm(p2vm)
        hdr = self.hdus[0][prod_name]["PRIMARY"].header
        hdr_ft = self.hdus[0][f"{prod_name}_raw_1"][ext_data_ft].header
        self.fname = self.hdus[0][prod_name].filename()
        self.procatg = hdr.get("HIERARCH ESO PRO CATG")

        p2 = Panel(n, 2, height_ratios=[1, 4])

        vspace = 0.2
        detector_list = list(self.hdus[0][f"{prod_name}_raw_1"][ext_detector_ft].data)
        dname_list = []
        dpos_list = []
        for i in range(len(detector_list)):
            dname = detector_list[i][4]
            dpos = detector_list[i][5][1]
            dname_list.append(dname)
            dpos_list.append(dpos)

        t1a = TextPlot(columns=1, v_space=vspace)
        col1a = (
            str(hdr.get("INSTRUME")),
            "EXTNAME: " + ext_data_ft, 
            "PRO CATG: "
            + str(hdr.get("HIERARCH ESO PRO CATG")),
            "FILE NAME: " + os.path.basename(self.fname),
            "RAW1.NAME: " + str(hdr.get("HIERARCH ESO PRO REC1 RAW1 NAME")),
            "NAXIS3: " + str(hdr_ft.get("NAXIS3")),
        )
        t1a.add_data(col1a)
        p2.assign_plot(t1a, 0, 0, xext=2)

        t2 = TextPlot(columns=1, v_space=vspace, xext=1)
        col2 = self.metadata[:4]
        t2.add_data(col2)
        p2.assign_plot(t2, 2, 0, xext=2)

        t3 = TextPlot(columns=1, v_space=vspace, xext=1)
        col3 = self.metadata[4:]
        t3.add_data(col3)
        p2.assign_plot(t3, 4, 0)
        
        if prod_name == "p2vm":
            nfiles = 6
        else:
            nfiles = 4

        #Raw image plots
        for i in range(nfiles):
            ft_data = self.hdus[0][f"{prod_name}_raw_{i+1}"][ext_data_ft].data
            arrays = [arr for _, arr in ft_data]
            im_data = np.stack(arrays, axis=0)
            a = np.shape(im_data.mean(axis=0).T)[1] / np.shape(im_data.mean(axis=0).T)[0]
            #Generate image plot
            im_plot= ImagePlot(
                im_data.mean(axis=0).T, 
                title=f"{prod_name.upper()}_RAW_{str(i+1)}",
                v_clip="percentile",
                v_clip_kwargs={"percentile": 95.0},
                #show_colorbar=False,
                aspect=2*a
            )
            p2.assign_plot(im_plot, i, 1)


        if include_product == 1:

            ft_data = self.hdus[0][prod_name][ext_data_ft].data[0][1]
            nx = np.shape(ft_data)[0]
            ny = np.shape(ft_data)[1]
            
            im_plot2 = ImagePlot(
                ft_data.T, 
                title=f"{self.procatg}",
                v_clip="percentile",
                v_clip_kwargs={"percentile": 95.0},
                #show_colorbar=False, 
                y_ticks=dpos_list,
                aspect = 2*(nx / ny)
                )
    
            p2.assign_plot(im_plot2, n-1, 1)

        return p2
        

    def panel_line_plots(self, sc_ft):
        """Function to generate line plots, 
        used to generate panels 7 to 12.

        Inputs: 
        sc_ft: Extension name of the raw data
        
        Outputs:
        Panel object p2 (image plots)
        """

        vspace = 0.2

        p2vm = self.hdus[0]["p2vm"]
        hdr = p2vm["PRIMARY"].header
        self.metadata = GravitySetupInfo.p2vm(p2vm)
        ext = [p2vm[i].header.get("EXTNAME") for i in range(1, len(p2vm))]
        #Wavelength data index for SC and FT
        index_sc = [idx for idx, s in enumerate(ext) if "OI_WAVELENGTH" in s][0]
        index_ft = [idx for idx, s in enumerate(ext) if "OI_WAVELENGTH" in s][-1]

        reg = p2vm[sc_ft].data["REGNAME"]   #Regname labels

        #Wavelength data and indexes for phase plots
        if sc_ft == "P2VM_SC":
            wave = [list(p2vm[index_sc+1].data)[i][0]/1e-6 for i in range(len(list(p2vm[index_sc+1].data)))]
            idx = [ 4, 2, 0, 5, 3, 1 ]
        else:
            wave = [list(p2vm[index_ft+1].data)[i][0]/1e-6 for i in range(len(list(p2vm[index_ft+1].data)))]
            idx = [ 5, 3, 4, 1, 2, 0 ]

        #Transmission, coherence, and phase data
        trans = p2vm[sc_ft].data['TRANSMISSION']
        coh = p2vm[sc_ft].data['COHERENCE']
        pi = 3.141593
        phase = p2vm[sc_ft].data['PHASE'] / 2.0 / pi * 360.0        
        if len(reg) == 48:
            nrows = 6
            mult = np.repeat(idx, 8)
        else:
            nrows = 3
            mult = np.repeat(idx, 4)

        #Metadata
        t1a = TextPlot(columns=1, v_space=vspace)
        col1a = (
            str(hdr.get("INSTRUME")),
            "EXTNAME: " + sc_ft, 
            "PRO CATG: "
            + str(hdr.get("HIERARCH ESO PRO CATG")),
            "FILE NAME: " + os.path.basename(self.fname),
            "RAW1.NAME: " + str(hdr.get("HIERARCH ESO PRO REC1 RAW1 NAME")),
            #"NAXIS3: " + str(hdr.get("NAXIS3")),
        )
        t1a.add_data(col1a)
        t2 = TextPlot(columns=1, v_space=vspace, xext=1)
        col2 = self.metadata[:4]
        t2.add_data(col2)

        t3 = TextPlot(columns=1, v_space=vspace, xext=1)
        col3 = self.metadata[4:]
        t3.add_data(col3)

        p5 = Panel(8, 7, height_ratios=[1, 4, 4, 4, 4, 4, 4])
        p5.assign_plot(t1a, 0, 0, xext=2)
        p5.assign_plot(t2, 2, 0, xext=1)
        p5.assign_plot(t3, 3, 0)

        p6 = Panel(8, 7, height_ratios=[1, 4, 4, 4, 4, 4, 4])
        p6.assign_plot(t1a, 0, 0, xext=2)
        p6.assign_plot(t2, 2, 0, xext=1)
        p6.assign_plot(t3, 3, 0)

        p7 = Panel(8, 7, height_ratios=[1, 4, 4, 4, 4, 4, 4])
        p7.assign_plot(t1a, 0, 0, xext=2)
        p7.assign_plot(t2, 2, 0, xext=1)
        p7.assign_plot(t3, 3, 0)

        color_arr = ["red", "blue", "green", "orange", "violet", "gray"]
        tel_arr = ["tel 1", "tel 2", "tel 3", "tel 4", "tel 5", "tel 6"]

        #Generate panels

        for k in range(nrows-1):

            for i in range(0, 8):
                tel_num = [int(reg[i+k*8][0])-1,int(reg[i+k*8][1])-1]
    
                trans_tel1 = trans[i+k*8,tel_num[0],:]
                trans_tel2 = trans[i+k*8,tel_num[1],:]
    
                total = trans[i+k*8,0,:]+trans[i+k*8,1,:]+trans[i+k*8,2,:]+trans[i+k*8,3,:]
                total2 = 2 * np.sqrt(abs(trans_tel1*trans_tel2))
                
                lp = LinePlot(legend = True,
                             y_label="rel. trans.",
                             title=reg[i+k*8])
                lp.set_ylim(ymin=0, ymax=1)
                for j in range(np.shape(trans)[1]):
                    lp.add_data(d=[wave, list(trans[i+k*8,j,:]/total)], color = color_arr[j], label=tel_arr[j])
    
                cp = LinePlot(legend = False,
                             y_label="norm. coher.",
                             title=reg[i+k*8])
                cp.set_ylim(ymin=0.7, ymax=1.1)
                for j in range(np.shape(coh)[1]):
                    cp.add_data(d=[wave, list(coh[i+k*8,j,:]/total2)], color = color_arr[j], label=tel_arr[j])
                
                pp = LinePlot(legend = False,
                             y_label="phase",
                             title=reg[i+k*8])
                            
                cur_phase = list(phase[i+k*8,mult[i+k*8],:])
                pp.set_ylim(ymin=-20.0, ymax = +380.0)
                pp.add_data(d=[wave, cur_phase], color = color_arr[mult[i+k*8]], label=tel_arr[mult[i+k*8]])
        
    
                p5.assign_plot(lp, i, k+1)
                p6.assign_plot(cp, i, k+1)
                p7.assign_plot(pp, i, k+1)
            

        for i in range(0, 8):
            tel_num = [int(reg[i+8*(nrows-1)][0])-1,int(reg[i+8*(nrows-1)][1])-1]

            trans_tel1 = trans[i+8*(nrows-1),tel_num[0],:]
            trans_tel2 = trans[i+8*(nrows-1),tel_num[1],:]


            total = trans[i+8*(nrows-1),0,:]+trans[i+8*(nrows-1),1,:]+trans[i+8*(nrows-1),2,:]+trans[i+8*(nrows-1),3,:]
            total2 = 2 * np.sqrt(abs(trans_tel1*trans_tel2))
            
            lp = LinePlot(legend = True,
                        x_label="lambda / mu",
                         y_label="rel. trans.",
                         title=reg[i+8*(nrows-1)])
            lp.set_ylim(ymin=0, ymax=1)
            for j in range(np.shape(trans)[1]):
                lp.add_data(d=[wave, list(trans[i+8*(nrows-1),j,:]/total)], color = color_arr[j], label=tel_arr[j])

            cp = LinePlot(legend = False,
                        x_label="lambda / mu",
                         y_label="norm. coher.",
                         title=reg[i+8*(nrows-1)])
            cp.set_ylim(ymin=0.7, ymax=1.1)

            for j in range(np.shape(coh)[1]):
                cp.add_data(d=[wave, list(coh[i+8*(nrows-1),j,:]/total2)], color = color_arr[j], label=tel_arr[j])

            pp = LinePlot(legend = False,
                            x_label="lambda / mu",
                         y_label="phase",
                         title=reg[i+8*(nrows-1)])
                        
            cur_phase = list(phase[i+8*(nrows-1),mult[i+8*(nrows-1)],:])
            pp.set_ylim(ymin=-20.0, ymax = +380.0)
            pp.add_data(d=[wave, cur_phase], color = color_arr[mult[i+8*(nrows-1)]], label=tel_arr[mult[i+8*(nrows-1)]])

            p5.assign_plot(lp, i, nrows)
            p6.assign_plot(cp, i, nrows)
            p7.assign_plot(pp, i, nrows)


        return p5, p6, p7


    
    def generate_panels(self, **kwargs):
        panels = {}

        p1, p2 = self.panel_images_cuts_sc(5, "flat", 1)
        addme1 = {
            "panel": "first panel",
            "report_name": f"gravity_FLAT_SC_images",
            "report_description": f"GRAVITY p2vm panel 1",
            "report_tags": [],
        }

        addme2 = {
            "panel": "second panel",
            "report_name": f"gravity_FLAT_SC_plots",
            "report_description": f"GRAVITY p2vm panel 2",
            "report_tags": [],
        }
        
        p3, p4 = self.panel_images_cuts_sc(6, "p2vm", 0)
        addme3 = {
            "panel": "third panel",
            "report_name": f"gravity_P2VM_SC_images",
            "report_description": f"GRAVITY p2vm panel 3",
            "report_tags": [],
        }

        addme4 = {
            "panel": "fourth panel",
            "report_name": f"gravity_P2VM_SC_plots",
            "report_description": f"GRAVITY p2vm panel 4",
            "report_tags": [],
        }

        p5 = self.panel_images_cuts_ft(5, "flat", 1)
        addme5 = {
            "panel": "fifth panel",
            "report_name": f"gravity_FLAT_FT_images",
            "report_description": f"GRAVITY p2vm panel 5",
            "report_tags": [],
        }

        p6 = self.panel_images_cuts_ft(6, "p2vm", 0)
        addme6 = {
            "panel": "sixth panel",
            "report_name": f"gravity_P2VM_FT_images",
            "report_description": f"GRAVITY p2vm panel 6",
            "report_tags": [],
        }

        p7, p9, p11 = self.panel_line_plots("P2VM_SC")
        addme7 = {
            "panel": "seventh panel",
            "report_name": f"gravity_P2VM_SC_plots_2",
            "report_description": f"GRAVITY p2vm panel 7",
            "report_tags": [],
        }
        addme9 = {
            "panel": "ninth panel",
            "report_name": f"gravity_P2VM_SC_plots_3",
            "report_description": f"GRAVITY p2vm panel 9",
            "report_tags": [],
        }
        addme11= {
            "panel": "eleventh panel",
            "report_name": f"gravity_P2VM_SC_plots_4",
            "report_description": f"GRAVITY p2vm panel 11",
            "report_tags": [],
        }        

        p8, p10, p12 = self.panel_line_plots("P2VM_FT")
        addme8 = {
            "panel": "eight panel",
            "report_name": f"gravity_P2VM_FT_plots_2",
            "report_description": f"GRAVITY p2vm panel 8",
            "report_tags": [],
        }
        addme10 = {
            "panel": "tenth panel",
            "report_name": f"gravity_P2VM_FT_plots_3",
            "report_description": f"GRAVITY p2vm panel 10",
            "report_tags": [],
        }
        addme12 = {
            "panel": "twelfth panel",
            "report_name": f"gravity_P2VM_FT_plots_4",
            "report_description": f"GRAVITY p2vm panel 12",
            "report_tags": [],
        }

    
        panels[p1] = addme1
        panels[p2] = addme2
        panels[p3] = addme3
        panels[p4] = addme4
        panels[p5] = addme5
        panels[p6] = addme6
        panels[p7] = addme7
        panels[p8] = addme8
        panels[p9] = addme9
        panels[p10] = addme10
        panels[p11] = addme11
        panels[p12] = addme12
            
        return panels


rep = GravityP2vmReport()
