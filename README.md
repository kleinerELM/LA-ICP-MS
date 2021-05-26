# 2D LA-ICP-MS
Script to process 2D LA-ICP-MS data from TUM.

It can process data structured as xlsx or can be provided as raw data. Examples can be found in the subfolder `example_data`.
This package provides a python script and a jupyter file. While the jupyter file should be self-explanatory, the usage of the python file is explained in the help dialog:

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

# Jupyter notebook

It is recommendet, to use the included jupyter notebook `laser.ipynb` for the analysis.

# required thrid party packages

This script was tested using Python 3.9.5

This package requires some third-party-packages, which can be installed via pip using this command:

```
pip install -r requirements.txt
```