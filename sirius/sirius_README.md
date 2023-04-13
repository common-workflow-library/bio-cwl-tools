# How to use SIRIUS in MAW:
Pass your SIRIUS user email and password and the ```.ms``` input file to the cwltool.
```
cwltool --enable-ext --cachedir cache --debug maw-sirius.cwl --sirius_user your_email --sirius_password your_password --spectrum /Users/mahnoorzulfiqar/Downloads/New_ms2_spectra_endo_pos/insilico/SIRIUS/no_isotope/1_NA_iso_NA_MS1p_446.230465650658_SIRIUS_param.ms
```