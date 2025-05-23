class GravityReportMixin(object):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._version = "TBD"


from . import gravity_master_dark
from . import gravity_wavelength
from . import gravity_rawdisp
from . import gravity_calibrator
from . import gravity_p2vm


__all__ = [
    gravity_master_dark,
    gravity_wavelength,
    gravity_rawdisp,
    gravity_calibrator,
    gravity_p2vm,
]
