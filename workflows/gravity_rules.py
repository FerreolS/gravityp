from edps import match

from . import gravity_keywords as kwd


def is_gravity(f):
    return f[kwd.instrume] == "GRAVITY"


def is_calib(f):
    return is_gravity(f) and f[kwd.dpr_catg] == "CALIB"


def is_science(f):
    return is_gravity(f) and f[kwd.dpr_catg] == "SCIENCE"


def is_single_sci(f):
    return is_science(f) and "OBJECT" in f[kwd.dpr_type] and "SINGLE" in f[kwd.dpr_type]


def is_single_science_sky(f):
    c1 = is_calib(f) and f[kwd.dpr_type] == "SKY,SINGLE" and f[kwd.tpl_id] != "GRAVITY_gen_tec_checkMetZero"
    c2 = is_science(f) and f[kwd.dpr_type] == "SKY,SINGLE"
    return c1 or c2


def is_dual_sci(f):
    return is_science(f) and "OBJECT" in f[kwd.dpr_type] and "DUAL" in f[kwd.dpr_type]


def is_dual_sky_calib(f):
    c1 = is_calib(f) and f[kwd.dpr_type] == "SKY,DUAL" and f[kwd.tpl_id] != "GRAVITY_gen_tec_checkMetZero"
    return c1


def is_single_sky_calib(f):
    c1 = is_calib(f) and f[kwd.dpr_type] == "SKY,SINGLE" and f[kwd.tpl_id] != "GRAVITY_gen_tec_checkMetZero"
    return c1


# ASSOCIATION RULES
#  -  first, e.g.  ref=trigger (e.g. science)
#  -  second, e.g.  f=file to associate (e.g. calibration)


def is_assoc(ref, f):
    return True


def is_assoc_dark_normal(ref, f):
    return (match(ref, f, [kwd.ins_polar_mode, kwd.ft_polar_mode, kwd.dit2, kwd.ins_spec_res, kwd.ins_fddl_window]) and
            f[kwd.tpl_id] == "GRAVITY_gen_cal_dark")


def is_assoc_dark_p2vm(ref, f):
    return (match(ref, f, [kwd.ins_polar_mode, kwd.ft_polar_mode, kwd.dit2, kwd.ins_spec_res, kwd.ins_fddl_window]) and
            f[kwd.tpl_id] == "GRAVITY_gen_cal_p2vm")
