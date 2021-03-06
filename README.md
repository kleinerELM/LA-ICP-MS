# 2D LA-ICP-MS

Script to process 2D LA-ICP-MS data from TUM.

It can process data structured as xlsx or can be provided as raw data. Examples can be found in the subfolder `example_data`.
This package provides a python script and a jupyter file. While the jupyter file should be self-explanatory, the usage of the python file is explained in the help dialog:

# direct script usage

The direct usage of the processing script is possible but not endorsed.

Help output:
```
#########################################################
# extract images from LA-ICP-MS data by TUM             #
#                                                       #
# © 2022 Florian Kleiner                                #
#   Bauhaus-Universität Weimar                          #
#   Finger-Institut für Baustoffkunde                   #
#                                                       #
#########################################################

usage: .\process_image.py [-h] [-o] [-r] [-x:] [-y:] [-d]
-h                : show this help
-o                : setting output directory name [processed]
-r                : load raw data instead of a xlsx file
-x                : de-noise the data in x direction, standard filter size: 2 px
-y                : change interpolation between the lines in x direction, standard: 6 lines

```

# Jupyter notebook

It is recommended, to use the included jupyter notebook `laser+edx.ipynb` for the analysis.

## define the scale

Use a SEM image acquired after the LA-ICP-MS ablation and measure the outer dimensions of the ablated area.

## combine EDX segmentation and LA-ICP-MS

To combine an EDX segmentation with the LA-ICP-MS Data, the binary image has to be transformed to the orientation and dimensions of the LA-ICP-MS dataset.
In the example data, the images were combined using SEM images before and after the LA-ICP-MS ablation.

## phase map alignment

There is an alignment exaplme in the `example_data` folder. This will be polished in the future.

# required thrid party packages

This script was tested using Python 3.9.5

This package requires some third-party-packages, which can be installed via pip using this command:

```
pip install -r requirements.txt
```

# citation

this script is published together with the paper _Combined LA-ICP-MS and SEM-EDX analyses for spatially resolved major minor and trace element detection in cement clinker phases, Cement and Concrete Research_

cite as following:

Florian Kleiner, Marco Decker, Christiane Rößler, Harald Hilbig, Horst-Michael Ludwig, _Combined LA-ICP-MS and SEM-EDX analyses for spatially resolved major minor and trace element detection in cement clinker phases, Cement and Concrete Research_, 2022, doi: [10.1016/j.cemconres.2022.106875](https://doi.org/10.1016/j.cemconres.2022.106875)