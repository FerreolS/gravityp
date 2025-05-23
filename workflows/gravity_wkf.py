from edps import qc1calib, qc0, calchecker, ReportInput
from edps import task

from .gravity_calibrate_visibilities import calibrate_visibilities
from .gravity_datasources import *
from .gravity_quality_control import quality_control

__title__ = "GRAVITY workflow"
idp = "idp"

# --- Processing tasks -----------------------------------------------------------------------------

# - Task that processes the raw darks.
calib_dark = (task('dark')
              .with_recipe('gravity_dark')
              .with_report('gravity_rawdisp', ReportInput.RECIPE_INPUTS)
              .with_report('gravity_master_dark', ReportInput.RECIPE_INPUTS_OUTPUTS)
              .with_main_input(raw_dark)
              .with_meta_targets([qc1calib])
              .build())

# - Task that processes pixel to visibility (P2VM) internal calibrations. It calibrates the instrument bad pixels,
#   wavelength table, interferometric contrast and phase.
pixel_to_visibility = (task('pixel_to_visibility')
                       .with_recipe('gravity_p2vm')
                       .with_report('gravity_rawdisp', ReportInput.RECIPE_INPUTS)
                       .with_report('gravity_p2vm', ReportInput.RECIPE_INPUTS_OUTPUTS)
                       .with_main_input(raw_calibrations)
                       .with_associated_input(calib_dark, [DARK_P2VM], match_rules=match_dark_p2vm)
                       .with_associated_input(wave_param)
                       .with_meta_targets([qc1calib, calchecker])
                       .build())

# - Subworkflow that runs tasks dedicated to quality control and instrument monitoring
actuators_response, arc, dispersion = quality_control(calib_dark, pixel_to_visibility)

# - Task that computes the uncalibrated visibilites for the calibrator
calibrator_visibilities = (task('calibrator_visibilities')
                           .with_recipe('gravity_vis')
                           .with_report('gravity_rawdisp', ReportInput.RECIPE_INPUTS)
                           .with_report('gravity_calibrator', ReportInput.RECIPE_INPUTS_OUTPUTS)
                           .with_main_input(raw_cal_object)
                           .with_associated_input(calib_dark, [DARK], match_rules=match_dark_normal)
                           .with_associated_input(pixel_to_visibility,
                                                  [MASTER_FLAT, MASTER_WAVE, MASTER_P2VM, MASTER_BAD])
                           .with_associated_input(disp_model)
                           .with_associated_input(static_param)
                           .with_associated_input(eop_param, min_ret=0)
                           .with_associated_input(diameter_catalog, min_ret=0)
                           .with_meta_targets([qc1calib, calchecker])
                           .build())

# - Task that computes the uncalibrated visibilites for the science target
target_visiblities = (task('target_visibilities')
                      .with_recipe('gravity_vis')
                      .with_main_input(raw_sci_object)
                      .with_associated_input(calib_dark, [DARK], match_rules=match_dark_normal)
                      .with_associated_input(pixel_to_visibility,
                                             [MASTER_FLAT, MASTER_WAVE, MASTER_P2VM, MASTER_BAD])
                      .with_associated_input(disp_model)
                      .with_associated_input(static_param)
                      .with_associated_input(eop_param, min_ret=0)
                      .with_associated_input(diameter_catalog, min_ret=0)
                      .with_meta_targets([qc0, calchecker])
                      .build())

# Subworkflow that computes the calibrated visibilities, either starting from the raw data or from pre-processed data
# (e.g., downloaded from the archive or on disk)
calibrate_visibilites, calibrate_pre_computed_visibilites = calibrate_visibilities(target_visiblities,
                                                                                   calibrator_visibilities)
