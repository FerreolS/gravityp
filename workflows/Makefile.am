AUTOMAKE_OPTIONS = foreign

#Put here all the workflows and associate python files here
WORKFLOWFILES = gravity_calibrate_visibilities.py gravity_classification.py gravity_datasources.py gravity_keywords.py gravity_quality_control.py gravity_rules.py gravity_wkf.py __init__.py

#Put here all the recipe parameters associated with the workflow
PARAMETERFILES = gravity_parameters.yaml

#This installs the workflow and the OCA rules in reflexaux directory
workflow_DATA = $(WORKFLOWFILES) ${PARAMETERFILES}

EXTRA_DIST = $(WORKFLOWFILES) ${PARAMETERFILES}
