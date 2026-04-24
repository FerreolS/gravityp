from edps import data_source, match_rules
from edps.generator.time_range import *

from .gravity_classification import *

# Data sources
# Convention for Data sources Association rule levels:
# Each data source can have several match function which correspond to different
# quality levels for the selected data. The level is specified as a number that
# follows this convention:
#   level < 0: more restrictive than the calibration plan
#   level = 0 follows the calibration plan
#   level = 1 quality sufficient for QC1 certification
#   level = 2 probably still acceptable quality
#   level = 3 significant risk of bad quality results

dit_keywords = [kwd.dit1, kwd.dit2, kwd.dit3]
polar_keywords = [kwd.ins_polar_mode, kwd.ft_polar_mode]
pointing_keywords = [kwd.ins_sobj_offx, kwd.ins_sobj_offy]
group_dark = dit_keywords + polar_keywords + [kwd.ins_spec_res, kwd.tpl_start, kwd.ins_fddl_window + kwd.tpl_id]
setup_dark = dit_keywords + polar_keywords + [kwd.ins_spec_res, kwd.ins_fddl_window]
match_dark = polar_keywords + [kwd.dit2, kwd.ins_spec_res, kwd.ins_fddl_window]
group_wave_lamp = dit_keywords + polar_keywords + [kwd.ins_spec_res, kwd.arcfile]
setup_wave_lamp = dit_keywords + polar_keywords + [kwd.ins_spec_res]
telescope_keywords = [kwd.iss_conf_station1, kwd.iss_conf_station2, kwd.iss_conf_station3, kwd.iss_conf_station4]
match_p2vm = polar_keywords + [kwd.ins_spec_res]
setup = polar_keywords + [kwd.dit2, kwd.dit3, kwd.ins_spec_res, kwd.gain]
group = polar_keywords + [kwd.dit2, kwd.dit3, kwd.ins_spec_res, kwd.gain, kwd.tpl_start]
group_science = polar_keywords + pointing_keywords + [kwd.dit2, kwd.dit3, kwd.ins_spec_res, kwd.gain, kwd.tpl_start]
setup_dispersion = polar_keywords + [kwd.dit1, kwd.dit2, kwd.dit3, kwd.ins_spec_res]

ONE_AND_HALF_DAYS = RelativeTimeRange(-1.5, 1.5)

# Raw types

# --- Datasources and association rules for dark ---
# There are 2 types of darks, they differ by the tpl.id header keywords
#   - "Normal" darks, with tpl.id = "GRAVITY_gen_cal_dark".
#   - "P2VM" darks, with tpl.id = "GRAVITY_gen_cal_p2vm" (P2VM stands for pixel to visibility map).
# They both are processed by the same recipe.
# Some tasks require normal darks, some task require p2vm darks, and other can take either one of them.
# For the first two case, match rules are defined below.
# Tasks without preference simply pick up the dark closest in time using the match rule defined in the data source.

raw_dark = (data_source("DARK")
            .with_classification_rule(dark_class)
            .with_setup_keywords(setup_dark)
            .with_grouping_keywords(group_dark)
            .with_match_keywords(match_dark, time_range=ONE_AND_HALF_DAYS, level=0)
            .with_match_keywords(match_dark, time_range=TWO_WEEKS, level=1)
            .with_match_keywords(match_dark, time_range=ONE_MONTH, level=3)
            .build())

match_dark_normal = (match_rules()
                     .with_match_function(rules.is_assoc_dark_normal, time_range=ONE_AND_HALF_DAYS, level=0)
                     .with_match_function(rules.is_assoc_dark_normal, time_range=TWO_WEEKS, level=1)
                     .with_match_function(rules.is_assoc_dark_normal, time_range=ONE_MONTH, level=3))

match_dark_p2vm = (match_rules()
                   .with_match_function(rules.is_assoc_dark_p2vm, time_range=ONE_AND_HALF_DAYS, level=0)
                   .with_match_function(rules.is_assoc_dark_p2vm, time_range=TWO_WEEKS, level=1)
                   .with_match_function(rules.is_assoc_dark_p2vm, time_range=ONE_MONTH, level=3))
# ----------------------------------------

raw_calibrations = (data_source("PIXEL_TO_VISIBILITY")
                    .with_classification_rule(flat_class)
                    .with_classification_rule(wave_class)
                    .with_classification_rule(wavesc_class)
                    .with_classification_rule(p2vm_class)
                    .with_setup_keywords(setup + [kwd.dit1])
                    .with_grouping_keywords(group + [kwd.dit1])
                    .with_match_keywords(match_p2vm, time_range=ONE_AND_HALF_DAYS, level=0)
                    .with_match_keywords(match_p2vm, time_range=TWO_WEEKS, level=1)
                    .with_match_keywords(match_p2vm, time_range=ONE_MONTH, level=3)
                    .build())

raw_wave_lamp = (data_source("ARC")  # check the name I used in other workflows.
                 .with_classification_rule(wave_lamp_class)
                 .with_grouping_keywords(group_wave_lamp)
                 .with_setup_keywords(setup_wave_lamp)
                 .with_match_keywords(match_p2vm, time_range=ONE_AND_HALF_DAYS, level=0)
                 .with_match_keywords(match_p2vm, time_range=TWO_WEEKS, level=1)
                 .with_match_keywords(match_p2vm, time_range=ONE_MONTH, level=3)
                 .build())

raw_piezotf = (data_source("PIEZO_ACTUACTORS")
               .with_classification_rule(piezotf_class)
               .with_grouping_keywords(group_wave_lamp)
               .with_setup_keywords(setup_wave_lamp)
               .with_match_keywords(match_p2vm, time_range=ONE_AND_HALF_DAYS, level=0)
               .with_match_keywords(match_p2vm, time_range=TWO_WEEKS, level=1)
               .with_match_keywords(match_p2vm, time_range=ONE_MONTH, level=3)
               .build())

raw_dispersion = (data_source("DISPERSION")
                  .with_classification_rule(dispersion_class)
                  .with_grouping_keywords(setup_dispersion + [kwd.tpl_start])
                  .with_setup_keywords(setup_dispersion)
                  .with_match_keywords(match_p2vm, time_range=ONE_AND_HALF_DAYS, level=0)
                  .with_match_keywords(match_p2vm, time_range=TWO_WEEKS, level=1)
                  .with_match_keywords(match_p2vm, time_range=ONE_MONTH, level=3)
                  .build())

# Standard calibrator
raw_cal_object = (data_source("CALIBRATOR")
                  .with_classification_rule(single_calib_class)
                  .with_classification_rule(dual_calib_class)
                  .with_setup_keywords(setup + telescope_keywords)
                  .with_grouping_keywords(group)
                  .with_match_keywords([kwd.ins_spec_res, kwd.ins_polar_mode, kwd.ft_polar_mode,
                                        kwd.telescop, kwd.night, "$associate_calibrators", "$match_exposure_time"],
                                       level=0, time_range=ONE_DAY)
                  .build())

raw_cal_sky = (data_source("CALIBRATOR_SKY")
               .with_classification_rule(single_sky_calib_class)
               .with_classification_rule(dual_sky_calib_class)
               .with_match_keywords([kwd.ins_spec_res, kwd.ins_polar_mode, kwd.ft_polar_mode,
                                     kwd.telescop, kwd.tpl_start], level=0, time_range=ONE_DAY)
               .build())

raw_sci_sky = (data_source("SKY")
               .with_classification_rule(single_sky_science_class)
               .with_classification_rule(dual_sky_science_class)
               .with_match_keywords([kwd.tpl_start])
               .build())

# Science target. Note that this is grouped by
# target within the same OB (equal INS.SOBJ.OFFX and INS.SOBJ.OFFY keywords)
# These keywords were added in 2017. Previous data with missing keywords
# will be processed with all the objects from the same TPL together.
raw_sci_object = (data_source("OBJECT")
                  .with_classification_rule(single_sci_class)
                  .with_classification_rule(dual_science_class)
                  .with_setup_keywords(setup + telescope_keywords)
                  .with_grouping_keywords(group_science)
                  .build())

# Pre-reduced uncalibrated visibility of the scientific target and calibrator (e.g. IDPs) downloaded from the archive.
target_pre_computed_visibilities = (data_source("TARGET_VISIBILITIES")
                                    .with_classification_rule(SINGLE_SCI_VIS)
                                    .with_classification_rule(DUAL_SCI_VIS)
                                    .with_setup_keywords(setup + telescope_keywords)
                                    .with_grouping_keywords([kwd.ins_spec_res, kwd.ins_polar_mode, kwd.pro_night_obs])
                                    .build())

calibrator_pre_computed_visibilites = (data_source("CALIBRATOR_VISIBILITIES")
                                       .with_classification_rule(SINGLE_CAL_VIS)
                                       .with_classification_rule(DUAL_CAL_VIS)
                                       .with_grouping_keywords(
    [kwd.ins_spec_res, kwd.ins_polar_mode, kwd.pro_night_obs])
                                       .with_match_keywords(
    [kwd.ins_spec_res, kwd.ins_polar_mode, kwd.telescop, kwd.ft_polar_mode, "night",
     "$associate_calibrators", "$match_exposure_time"])
                                       .build())

# Static calibrations
disp_model = (data_source()
              .with_classification_rule(disp_model_class)
              .with_match_keywords([kwd.instrume], time_range=IN_THE_PAST, level=0)
              .with_match_keywords([kwd.instrume], time_range=UNLIMITED, level=3)
              .build())

static_param = (data_source()
                .with_classification_rule(static_param_class)
                .with_match_keywords([kwd.instrume], time_range=IN_THE_PAST, level=0)
                .with_match_keywords([kwd.instrume], time_range=UNLIMITED, level=3)
                .build())

wave_param = (data_source()
              .with_classification_rule(wave_param_class)
              .with_match_keywords([kwd.instrume], time_range=IN_THE_PAST, level=0)
              .with_match_keywords([kwd.instrume], time_range=UNLIMITED, level=3)
              .build())

eop_param = (data_source()
             .with_classification_rule(eop_param_class)
             .with_match_function(rules.is_assoc).build())

diameter_catalog = (data_source("CALIBRATOR_CATALOG")
                    .with_classification_rule(diameter_catalogue_class)
                    .with_match_function(rules.is_assoc).build())
