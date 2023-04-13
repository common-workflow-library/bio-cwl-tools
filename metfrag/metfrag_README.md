# How to use MetFrag in MAW:
Pass the peaklist path, precursor mass, ion mode (1 for positive and 0 for negative), provide path to your local database, and a name for the result file
```
cwltool metfrag.cwl --PeakList path/to/peaklist.txt --IonizedPrecursorMass 123.45 --PrecursorIonMode 1 --IsPositiveIonMode True --LocalDatabase path/to/localdatabase.csv --SampleName result_file_name
```
