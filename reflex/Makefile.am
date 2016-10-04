AUTOMAKE_OPTIONS = foreign

#Put here all the workflows. They must end with extension .in. They should
#be like a normal workflow but using @prefix@ for paths 
#to the OCA rules (grep for OCA in the xml file) or Python scripts
WORKFLOWS = gravity_viscal.xml gravity_wkf.xml gravity_mastercalib.xml
#Put here the OCA rules used by the workflow
OCAWKF =  gravity_viscal_wkf.oca gravity_wkf.oca gravity_mastercalib.oca
#Put Python scripts used by the workflows
PYTHONWKF = gravity_plot_common.py gravity_viscal_interact.py

#This installs the workflow and the OCA rules in reflexaux directory
wkfextra_DATA = $(WORKFLOWS) $(OCAWKF) $(PYTHONWKF)

#Pack in the distribution the original workflow
EXTRA_DIST = $(WORKFLOWS).in $(OCAWKF) $(PYTHONWKF)

#This also installs the workflow in the common reflex directory
wkfcopy_DATA = $(WORKFLOWS)

