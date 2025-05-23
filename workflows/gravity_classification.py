from edps import classification_rule

from . import gravity_keywords as kwd
from . import gravity_rules as rules

gravity = {kwd.instrume: "GRAVITY"}
calibrations_kwd = {**gravity, kwd.pro_catg: None, kwd.dpr_catg: "CALIB"}
science_kwd = {**gravity, kwd.pro_catg: None, kwd.dpr_catg: "SCIENCE"}

# Raw types

dark_class = classification_rule('DARK_RAW', {**calibrations_kwd, kwd.dpr_type: "DARK"})
flat_class = classification_rule('FLAT_RAW', {**calibrations_kwd, kwd.dpr_type: "FLAT"})
wave_class = classification_rule('WAVE_RAW', {**calibrations_kwd, kwd.dpr_type: "WAVE"})

wavesc_class = classification_rule('WAVESC_RAW',
                                   {**calibrations_kwd,
                                    kwd.dpr_tech: ["INTERFEROMETRY,DIRECT", "INTERFEROMETRY,WOLLASTON"],
                                    kwd.dpr_type: ["WAVE,SC", "WAVESC"], kwd.tpl_id: "GRAVITY_gen_cal_p2vm"})

wave_lamp_class = classification_rule('WAVELAMP_RAW',
                                      {**calibrations_kwd, kwd.dpr_tech: ["SPECTRUM,DIRECT", "SPECTRUM,WOLLASTON"],
                                       kwd.dpr_type: ["WAVE,LAMP", "WAVELAMP"]})

piezotf_class = classification_rule('PIEZOTF_RAW', {**calibrations_kwd, kwd.dpr_type: "PIEZOTF"})
dispersion_class = classification_rule('DISP_RAW', {**calibrations_kwd, kwd.dpr_type: "DISP"})
p2vm_class = classification_rule('P2VM_RAW', {**calibrations_kwd, kwd.dpr_type: "P2VM"})
single_sci_class = classification_rule('SINGLE_SCI_RAW', rules.is_single_sci)
single_sky_science_class = classification_rule('SINGLE_SKY_RAW', {**science_kwd, kwd.dpr_type: "SKY,SINGLE"})
dual_science_class = classification_rule('DUAL_SCI_RAW', rules.is_dual_sci)
dual_sky_science_class = classification_rule('DUAL_SKY_RAW', {**science_kwd, kwd.dpr_type: "SKY,DUAL"})

single_calib_class = classification_rule('SINGLE_CAL_RAW',
                                         {**calibrations_kwd, kwd.dpr_type: ["OBJECT,SINGLE", "STD,SINGLE"]})
single_sky_calib_class = classification_rule('SINGLE_SKY_RAW', rules.is_single_sky_calib)
dual_calib_class = classification_rule('DUAL_CAL_RAW', {**calibrations_kwd, kwd.dpr_type: ["OBJECT,DUAL", "STD,DUAL"]})
dual_sky_calib_class = classification_rule('DUAL_SKY_RAW', rules.is_dual_sky_calib)

# Static calibrations
disp_model_class = classification_rule('DISP_MODEL', {**gravity, kwd.pro_catg: "DISP_MODEL"})
static_param_class = classification_rule('STATIC_PARAM', {**gravity, kwd.pro_catg: "STATIC_PARAM"})
wave_param_class = classification_rule('WAVE_PARAM', {**gravity, kwd.pro_catg: "WAVE_PARAM"})
eop_param_class = classification_rule("EOP_PARAM", {kwd.pro_catg: "EOP_PARAM"})
diameter_catalogue_class = classification_rule("DIAMETER_CAT", {kwd.pro_catg: "DIAMETER_CAT"})
# Master calibrations
DARK = classification_rule('DARK', {**gravity, kwd.pro_catg: 'DARK', kwd.tpl_id: "GRAVITY_gen_cal_dark"})
DARK_P2VM = classification_rule('DARK', {**gravity, kwd.pro_catg: 'DARK', kwd.tpl_id: "GRAVITY_gen_cal_p2vm"})
DARK_ANY = classification_rule('DARK', {**gravity, kwd.pro_catg: 'DARK'})

MASTER_FLAT = classification_rule('FLAT', {**gravity, kwd.pro_catg: 'FLAT'})
MASTER_WAVE = classification_rule('WAVE', {**gravity, kwd.pro_catg: 'WAVE'})
MASTER_P2VM = classification_rule('P2VM', {**gravity, kwd.pro_catg: 'P2VM'})
WAVELAMP_WKF = classification_rule('WAVELAMP', {**gravity, kwd.pro_catg: 'WAVELAMP'})
MASTER_BAD = classification_rule('BAD', {**gravity, kwd.pro_catg: 'BAD'})
DIODE_POSITION = classification_rule("DIODE_POSITION", {**gravity, kwd.pro_catg: "DIODE_POSITION"})

# Pre-reduced visibilities downloaded from archive (e.g. IDP)
SINGLE_CAL_VIS = classification_rule("SINGLE_CAL_VIS", {**gravity, kwd.pro_catg: "SINGLE_CAL_VIS"})
SINGLE_SCI_VIS = classification_rule("SINGLE_SCI_VIS", {**gravity, kwd.pro_catg: "SINGLE_SCI_VIS"})
DUAL_CAL_VIS = classification_rule("DUAL_CAL_VIS", {**gravity, kwd.pro_catg: "DUAL_CAL_VIS"})
DUAL_SCI_VIS = classification_rule("DUAL_SCI_VIS", {**gravity, kwd.pro_catg: "DUAL_SCI_VIS"})
