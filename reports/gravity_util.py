from adari_core.utils.utils import fetch_kw_or_default
import astropy.io.fits as fits

class GravityReportMixin(object):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._version = "TBD"

extensions = ["IMAGING_DATA_SC", "IMAGING_DATA_FT"]

def data_reader(filename):
    hdu = fits.open(filename, mode="readonly")
    ext1 = extensions[0]
    if len(hdu[ext1].data.shape) > 2:
        hdu[ext1].data = hdu[ext1].data[0]
    ext1 = extensions[1]
    hdu[ext1] = fits.ImageHDU(
        data=hdu[ext1].data[0][1].T, header=hdu[ext1].header, name=ext1
    )
    return hdu


def meta(hdul):
    metadata = [
        "INS.SPEC.RES: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO INS SPEC RES", default="N/A"
            )
        ),
        "INS.POLA.MODE: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO INS POLA MODE", default="N/A"
            )
        ),
        "FT.POLA.MODE: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO FT POLA MODE", default="N/A"
            )
        ),
        "DET1.SEQ1.DIT: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO DET1 SEQ1 DIT", default="N/A"
            )
        ),
        "DET2.SEQ1.DIT: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO DET2 SEQ1 DIT", default="N/A"
            )
        ),
        "DET3.SEQ1.DIT: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO DET3 SEQ1 DIT", default="N/A"
            )
        ),
        "INS.DET3.GAIN: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO INS DET3 GAIN", default="N/A"
            )
        ),
        "INS.MET.MODE: "
        + str(
            fetch_kw_or_default(
                hdul["PRIMARY"], "HIERARCH ESO INS MET MODE", default="N/A"
            )
        ),
    ]
    return metadata


class GravitySetupInfo:
    @staticmethod
    def dark(hdul):
        metadata = meta(hdul)
        return metadata

    @staticmethod
    def calib(hdul):
        metadata = meta(hdul)
        metadata = metadata + [
            "ISS.CONF.STATION1: "
            + str(
                fetch_kw_or_default(
                    hdul["PRIMARY"], "HIERARCH ESO ISS CONF STATION1", default="N/A"
                )
            ),
            "ISS.CONF.STATION2: "
            + str(
                fetch_kw_or_default(
                    hdul["PRIMARY"], "HIERARCH ESO ISS CONF STATION2", default="N/A"
                )
            ),
            "ISS.CONF.STATION3: "
            + str(
                fetch_kw_or_default(
                    hdul["PRIMARY"], "HIERARCH ESO ISS CONF STATION3", default="N/A"
                )
            ),
            "ISS.CONF.STATION4: "
            + str(
                fetch_kw_or_default(
                    hdul["PRIMARY"], "HIERARCH ESO ISS CONF STATION4", default="N/A"
                )
            ),
        ]
        return metadata

    @staticmethod
    def pixel_to_visibility(hdul):
        metadata = meta(hdul)
        return metadata

    @staticmethod
    def calibrator_visibilities(hdul):
        metadata = meta(hdul)
        return metadata

    @staticmethod
    def arc(hdul):
        metadata = meta(hdul)
        metadata = metadata[:-2]
        return metadata

    @staticmethod
    def actuators_response(hdul):
        metadata = meta(hdul)
        metadata = metadata[:-2]
        return metadata

    @staticmethod
    def p2vm(hdul):
        metadata = meta(hdul)
        return metadata
