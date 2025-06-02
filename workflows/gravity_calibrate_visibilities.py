from edps import subworkflow, science, calchecker, task

from .gravity_datasources import *

idp = "idp"


@subworkflow("calibrate_visibilities", "")
def calibrate_visibilities(target_visiblities, calibrator_visibilities):
    # This subworkflow computes the calibrated visibilities, either starting from raw data of target and calibrator
    # (task "calibrate_visibilities"), or pre-processed data for target and calibrator that are already present on
    # disk, e.g. downloaded from the ESO Science Portal or generated in a previous reduction (task
    # "calibrate_pre-computed_visibilities").

    # - This task creates calibrated visibilites out of the reduction of raw data.
    calibrate_visibilites = (task('calibrate_visibilities')
                             .with_recipe('gravity_viscal')
                             .with_main_input(target_visiblities)
                             .with_associated_input(calibrator_visibilities, min_ret=1, max_ret=100)
                             .with_associated_input(diameter_catalog, min_ret=0)
                             .with_meta_targets([idp, science, calchecker])
                             .build())

    # - This task creates calibrated visibilites from uncalibrated visibilites (e.g. downloaded from the
    # ESO science archive, or generated from previous reductions).
    calibrate_pre_computed_visibilites = (task('calibrate_pre_computed_visibilities')
                                          .with_recipe('gravity_viscal')
                                          .with_main_input(target_pre_computed_visibilities)
                                          .with_associated_input(calibrator_pre_computed_visibilites, min_ret=1,
                                                                 max_ret=20)
                                          .with_associated_input(diameter_catalog, min_ret=0)
                                          .with_meta_targets([idp, science])
                                          .build())

    return calibrate_visibilites, calibrate_pre_computed_visibilites
