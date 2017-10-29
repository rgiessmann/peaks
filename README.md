# peaks

## To execute run_* files:
* Download this repository and unzip it to `C:\peaks`
* Open "anaconda prompt"
* navigate to the folder containing the file called "input_traces.csv" (you can do this by typing: `cd ...` where ... can be the absolute path, e.g. `C:\Users\Robert\Desktop\`)
* type: `ipython C:\peaks\scripts\run_*.ipy` , where you replace * with the filename you want to run

## To generate annotation-able csv files
* Check https://gitlab.tu-berlin.de/r.giessmann/fsa_to_csv_converter

## To annotate .csv files of analyzed peaks with raw peak heights
* Open "anaconda prompt"
* navigate to the folder containing the file called "input_traces.csv" (you can do this by typing: `cd ...` where ... can be the absolute path, e.g. `C:\Users\Robert\Desktop\`)
* type: `annotate -i analyzed_peaks.csv ` , where you replace * with the filename you want to run



# Data specifications

## Input files

### input_traces.csv

tab separated values

|File name |Dye color | Ltot concentration (µM)
| --- | --- | ---
|01-18-16-11-27 AM.fsa|B|0
|01-18-16-35-11 AM.fsa|B|5

--> auto-picking of reference trace from all 0M traces?

### referenced .csv files of analyzed peaks

comma separated values

| Status|Sample Name|Sample Type|Size Standard|Analysis Method|SQI|Offscale|Quality|UD1|UD2|UD3|Dye|Sample Peak|Sample File Name|Size|Height|Area in Point|Area in BP|Data Point|Begin Point|Begin BP|End Point|End BP|Width in Point|Width in BP|User Comments|User Edit
| --- | --- | ---| --- | --- | ---| --- | --- | ---| --- | --- | ---| --- | --- | ---| --- | --- | ---| --- | --- | ---| --- | --- | ---| --- | --- | ---| --- | --- | ---
Analyzed|0|Sample|bto 550|Sizing Default - NPP|false|Pass|Pass||||B|1|01-18-16-11-27 AM.fsa||123|1117|0|2065|2052|-1|2072|-1|10|0||
Analyzed|0|Sample|bto 550|Sizing Default - NPP|false|Pass|Pass||||B|2|01-18-16-11-27 AM.fsa||83|911|0|2079|2072|-1|2091|-1|10|0||
