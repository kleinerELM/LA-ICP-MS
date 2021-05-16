# LA-ICP-MS
Script to process LA-ICP-MS data from TUM

It can process data structured as xlsx or can be provided as raw data. Examples can be found in the subfolder `example_data`


# direct script usage
```
#########################################################
# extract images from LA-ICP-MS data by TUM             #
#                                                       #
# © 2021 Florian Kleiner                                #
#   Bauhaus-Universität Weimar                          #
#   Finger-Institut für Baustoffkunde                   #
#                                                       #
#########################################################

usage: .\process_image.py [-h] [-o] [-r] [-x:] [-y:] [-d]
-h                : show this help
-o                : setting output directory name [processed]
-r                : load raw data instead of a xlsx file
-x                : denoise the data in x direction, standard filtersize: 2 px
-y                : change interpolation between the lines in x direction, standard: 6 lines

```
