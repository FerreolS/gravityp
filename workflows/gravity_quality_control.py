from edps import subworkflow, task, qc1calib, ReportInput

from .gravity_datasources import *


@subworkflow("quality_control", "")
def quality_control(calib_dark, pixel_to_visibility):
    # - Task to compute the response (open loop transfer function) of the piezo actuators used to fringe-track.
    actuators_response = (task('actuators_response')
                          .with_recipe('gravity_piezo')
                          .with_report('gravity_rawdisp', ReportInput.RECIPE_INPUTS)
                          .with_main_input(raw_piezotf)
                          .with_meta_targets([qc1calib])
                          .build())

    # - Task to reduce the raw Argon calibration lamp lines, which will be later used to compute the fiber dispersion
    arc = (task('arc')
           .with_recipe('gravity_wavelamp')
           .with_report('gravity_rawdisp', ReportInput.RECIPE_INPUTS)
           .with_report('gravity_wavelength', ReportInput.RECIPE_INPUTS_OUTPUTS)
           .with_main_input(raw_wave_lamp)
           .with_associated_input(calib_dark, [DARK_ANY])
           .with_associated_input(pixel_to_visibility, [MASTER_FLAT, MASTER_WAVE, MASTER_P2VM])
           .with_associated_input(wave_param)
           .with_meta_targets([qc1calib])
           .build())

    # - Task to measures the phases obtained on reduced Argon lines calibrations (task "arc" )
    #   and to measure the various positions (i.e., the fiber stretch) of the Fibered Differential Delay Lines (FDDL).
    #   The task computes the linearity model and the dispersion model of the differential delay lines.
    #   These models are used for quality control only, and not in the data reduction cascade.
    dispersion = (task('dispersion')
                  .with_recipe('gravity_disp')
                  .with_main_input(raw_dispersion)
                  .with_associated_input(calib_dark, [DARK_ANY])
                  .with_associated_input(pixel_to_visibility, [MASTER_FLAT, MASTER_WAVE, MASTER_P2VM, MASTER_BAD])
                  .with_associated_input(arc, [WAVELAMP_WKF])
                  .with_associated_input(static_param)
                  .with_meta_targets([qc1calib])
                  .build())

    return actuators_response, arc, dispersion
