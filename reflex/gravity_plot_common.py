try:
  import numpy
  import math
  from pipeline_display import *
  from pipeline_product import *
  from reflex_plot_widgets import *
  numpy.seterr(invalid='ignore')
except ImportError:
  donothing=1

class GravityInstrumentConf:
  def __init__(self, conf_fits):
    #This is for the time being hard-coded. It relates the gravity input
    #which runs from 1 to 4 witht he VLTI laboratory input. This will eventually
    #be part of the OPTICAL_TRAIN input.
    #This is also hard-coded in the pipeline (file gravi_utils.h, 
    #GRAVI_LABINPUT_* defines)
    self.lab_input = {}
    self.lab_input[1] =7
    self.lab_input[2] =5
    self.lab_input[3] =3
    self.lab_input[4] =1
    self.conf = PipelineProduct(conf_fits)

  def getLabInputIndex(self, gravity_input):
    return self.lab_input[gravity_input]
  
  def getStationName(self, gravity_input):
    this_lab_input = self.getLabInputIndex(gravity_input)
    ncol_ot = self.conf.getTableNcols('OPTICAL_TRAIN')
    tel_name_col = self.conf.readTableColumn('OPTICAL_TRAIN', 'TEL_NAME')
    for col in range(1, (ncol_ot - 2 ) / 2 + 1) :
      colname = "OPTI_NAME{0}".format(col)
      opti_elem = self.conf.readTableColumn('OPTICAL_TRAIN', colname)
      if(opti_elem[0] == 'LAB_INPUT1'):
        colname = "VALUE{0}".format(col)
        lab_input_col = self.conf.readTableColumn('OPTICAL_TRAIN', colname)
        row = numpy.where(lab_input_col == this_lab_input)[0][0]
        tel_name = tel_name_col[row]
        pass
    
    tel_col_geom = self.conf.readTableColumn('ARRAY_GEOMETRY', 'TEL_NAME')
    row_tel = numpy.where(tel_col_geom == tel_name)[0][0]
    sta_name_col = self.conf.readTableColumn('ARRAY_GEOMETRY', 'STA_NAME')
    sta_name = sta_name_col[row_tel]
    return sta_name

class PlotableRawCalibratorTF:
  def __init__(self, calib_list_fits):
    self.calib_list = list()
    for fits in calib_list_fits:
      self.calib_list.append(PipelineProduct(fits))
    self.loadFromFits()
    self.tfdisp = SpectrumDisplay()
    self.tfscdisp = ScatterDisplay()

  def loadFromFits(self):
    refhead = self.calib_list[0].all_hdu[0].header
    pola_mode = refhead['ESO INS POLA MODE']
    if pola_mode == 'SPLIT':
      self.npola = 2
    else:
      self.npola = 1
    self.gravi_input_sta = {}
    if 'ESO QC GRAVI_INPUT1 STA' in refhead:
      self.gravi_input_sta[1] = refhead['ESO QC GRAVI_INPUT1 STA']
    else:
      self.gravi_input_sta[1] = 'INP1'
    if 'ESO QC GRAVI_INPUT2 STA' in refhead:
      self.gravi_input_sta[2] = refhead['ESO QC GRAVI_INPUT2 STA']
    else:
      self.gravi_input_sta[2] = 'INP2'
    if 'ESO QC GRAVI_INPUT3 STA' in refhead:
      self.gravi_input_sta[3] = refhead['ESO QC GRAVI_INPUT3 STA']
    else:
      self.gravi_input_sta[3] = 'INP3'
    if 'ESO QC GRAVI_INPUT4 STA' in refhead:
      self.gravi_input_sta[4] = refhead['ESO QC GRAVI_INPUT4 STA']
    else:
      self.gravi_input_sta[4] = 'INP4'

    #The transfer functions to define:
    self.tf = {}
    self.time = list()
    #Loop on possible combinations of inputs and polarisations
    for pola in range(1, self.npola + 1):
      for input1 in range(1,4 + 1):
        for input2 in range(input1+1,4 + 1):
          confkey_sc = "SC{0}{1}_P{2}".format(input1, input2, pola)
          confkey_ft = "FT{0}{1}_P{2}".format(input1, input2, pola)
          self.tf[confkey_sc] = list()
          self.tf[confkey_ft] = list()
          #Loop on calibrators
          for calib in self.calib_list:
            head = calib.all_hdu[0].header
            keyname = "ESO QC TF VIS2_{0} MED".format(confkey_sc)
            value_sc = head[keyname]
            keyname = "ESO QC TF VIS2_{0} MED".format(confkey_ft)
            value_ft = head[keyname]
            self.tf[confkey_sc].append(value_sc) 
            self.tf[confkey_ft].append(value_ft) 

    #Get the timestamps
    for calib in self.calib_list:
      head = calib.all_hdu[0].header
      self.time.append(head['MJD-OBS']) 

    #Add two more points, at the beginning and the end, 
    #to avoid the degenerat case of one single calibrator
    mintime = min(self.time)
    maxtime = max(self.time)
    mintimeidx = numpy.argmin(self.time)
    maxtimeidx = numpy.argmax(self.time)
    self.time.insert(0, math.floor(mintime-0.5)+0.5)
    self.time.append(math.ceil(maxtime+.5)-.5)
    for pola in range(1, self.npola + 1):
      for input1 in range(1,4 + 1):
        for input2 in range(input1+1,4 + 1):
          confkey_sc = "SC{0}{1}_P{2}".format(input1, input2, pola)
          confkey_ft = "FT{0}{1}_P{2}".format(input1, input2, pola)
          mintfsc = self.tf[confkey_sc][mintimeidx]
          maxtfsc = self.tf[confkey_sc][maxtimeidx]
          self.tf[confkey_sc].insert(0, mintfsc)
          self.tf[confkey_sc].append(maxtfsc)
          mintfft = self.tf[confkey_ft][mintimeidx]
          maxtfft = self.tf[confkey_ft][maxtimeidx]
          self.tf[confkey_ft].insert(0, mintfft)
          self.tf[confkey_ft].append(maxtfft)

    #Set the colors
    self.color = {}
    self.color['SC12_P1'] = 'green'
    self.color['SC13_P1'] = 'blue'
    self.color['SC14_P1'] = 'red'
    self.color['SC23_P1'] = 'yellow'
    self.color['SC24_P1'] = 'orange'
    self.color['SC34_P1'] = 'magenta'
    self.color['SC12_P2'] = 'darkgreen'
    self.color['SC13_P2'] = 'darkblue'
    self.color['SC14_P2'] = 'darkred'
    self.color['SC23_P2'] = 'lightyellow'
    self.color['SC24_P2'] = 'darkorange'
    self.color['SC34_P2'] = 'darkmagenta'
    self.color['FT12_P1'] = 'green'
    self.color['FT13_P1'] = 'blue'
    self.color['FT14_P1'] = 'red'
    self.color['FT23_P1'] = 'yellow'
    self.color['FT24_P1'] = 'orange'
    self.color['FT34_P1'] = 'magenta'
    self.color['FT12_P2'] = 'darkgreen'
    self.color['FT13_P2'] = 'darkblue'
    self.color['FT14_P2'] = 'darkred'
    self.color['FT23_P2'] = 'lightyellow'
    self.color['FT24_P2'] = 'darkorange'
    self.color['FT34_P2'] = 'darkmagenta'

    
      
  def plot(self, subplots, title, tooltip):
    tf_1 = list()
    for i in range(len(self.time)):
      tf_1.append(1)
    isub = 0
    for input1 in range(1,4 + 1):
      
      for input2 in range(input1+1,4 + 1):
        if isub == 0 or isub == 1:
          this_title = title
        else:
          this_title = ''
        if isub == 4 or isub == 5:
          this_xlab = 'MJD-OBS'
        else:
          this_xlab = ''
        if isub == 2 :
          this_ylab = 'Transfer function of squared visibilities'
        else:
          this_ylab = ''
        self.tfdisp.setLabels(this_xlab, this_ylab)
        self.tfdisp.display(subplots[isub], this_title, tooltip, self.time, tf_1)
        labels = list()
        labels.append('vis=1')
        for pola in range(1, self.npola + 1):
          confkey_sc = "SC{0}{1}_P{2}".format(input1, input2, pola)
          confkey_ft = "FT{0}{1}_P{2}".format(input1, input2, pola)
          self.tfdisp.overplot(subplots[isub], self.time, 
                              self.tf[confkey_sc], self.color[confkey_sc]) 
          self.tfdisp.overplot(subplots[isub], self.time, 
                              self.tf[confkey_ft], self.color[confkey_ft]) 
          xlim1, xlim2 = self.tfdisp.wave_lim
          ylim1, ylim2 = self.tfdisp.flux_lim
          self.tfscdisp.setLimits(xlim1, xlim2, ylim1, ylim2)
          self.tfscdisp.setPointSize(10)
          self.tfscdisp.setColor(self.color[confkey_sc])
          self.tfscdisp.display(subplots[isub], '', '', self.time, self.tf[confkey_sc])
          self.tfscdisp.setColor(self.color[confkey_ft])
          self.tfscdisp.display(subplots[isub], '', '', self.time, self.tf[confkey_ft])

          if(self.npola) == 1:
            label_sc = "SC {0} - {1} ".format(self.gravi_input_sta[input1], 
                                              self.gravi_input_sta[input2])
            label_ft = "FT {0} - {1} ".format(self.gravi_input_sta[input1], 
                                              self.gravi_input_sta[input2])
          else:
            label_sc = "SC {0} - {1} Pol {2}".format(self.gravi_input_sta[input1], 
                                                     self.gravi_input_sta[input2],
                                                     pola)
            label_ft = "FT {0} - {1} Pol {2}".format(self.gravi_input_sta[input1], 
                                                     self.gravi_input_sta[input2],
                                                     pola)

          labels.append(label_sc)
          labels.append(label_ft)
          subplots[isub].legend(labels, loc='best', prop={'size':8})
        isub = isub + 1
