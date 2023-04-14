# How to use SIRIUS in MAW:
Pass your SIRIUS user email and password and the ```.ms``` input file to the cwltool.
```
cwltool --enable-ext --cachedir cache --debug maw-sirius.cwl --sirius_user your_email --sirius_password your_password --spectrum /path/to/SIRIUS_param.ms
```
