# -*- coding: utf-8 -*-

from PIL import Image
import os, sys, re, napari, getopt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import tkinter as tk
from tkinter import filedialog

def programInfo():
    print("#########################################################")
    print("# extract images from LA-ICP-MS data by TUM             #")
    print("#                                                       #")
    print("# © 2021 Florian Kleiner                                #")
    print("#   Bauhaus-Universität Weimar                          #")
    print("#   Finger-Institut für Baustoffkunde                   #")
    print("#                                                       #")
    print("#########################################################")
    print()

# Initial function to load the settings
def getBaseSettings():
    settings = {
        "showDebuggingOutput"  : False,
        "home_dir"             : os.path.dirname(os.path.realpath(__file__)),
        "workingDirectory"     : "",
        "excel_file"           : "",
        "outputDirectory"      : "processed",
        # gauss_smoothing_factor search a good value from 1 to 7 delivers good results. to disable smoothing, set the value to 0
        "smooth_y"             : 2,
        # interpolate lines between the lines in x direction (linear interpolation)
        "stretch_x"            : 6,
        "load_raw"             : False,
        # extracted from the xlsx - is it constant?
        "spot_distance_y"      : 0.579150579150579
    }
    return settings

#### process given command line arguments
def processArguments():
    settings = getBaseSettings()
    argv = sys.argv[1:]
    usage = sys.argv[0] + " [-h] [-o] [-r] [-x:] [-y:] [-d]"
    try:
        opts, args = getopt.getopt(argv,"horx:y:d",[])
    except getopt.GetoptError:
        print( usage )
    for opt, arg in opts:
        if opt == '-h':
            print( 'usage: ' + usage )
            print( '-h                : show this help' )
            print( '-o                : setting output directory name [{}]'.format(settings["outputDirectory"]) )
            print( '-r                : load raw data instead of a xlsx file' )
            print( '-x                : denoise the data in x direction, standard filtersize: {} px'.format(settings["smooth_y"]) )
            print( '-y                : change interpolation between the lines in x direction, standard: {} lines'.format(settings["stretch_x"]) )
            print( '' )
            sys.exit()
        elif opt in ("-o"):
            settings["outputDirectory"] = arg
            print( 'changed output directory to {}'.format(settings["outputDirectory"]) )
        elif opt in ("-r"):
            settings["load_raw"] = True
        elif opt in ("-y"):
            settings["smooth_y"] = int( arg )
            print( 'gaussian filter in y direction is set to {}'.format(settings["smooth_y"]) )
        elif opt in ("-x"):
            settings["stretch_x"] = int( arg )
            print( 'interpolate [{}] values between the lines in x direction (linear interpolation)'.format(settings["stretch_x"]) )
        elif opt in ("-d"):
            print( 'show debugging output' )
            settings["showDebuggingOutput"] = True
    print( '' )
    return settings


class LA_ICP_MS_LOADER:
    images      = {}
    np_images   = {}
    elements    = {}
    element_max = {} # max value in the data
    colormaps   = ["red" , "green", "blue", "yellow", "magenta", "cyan", "bop blue", "bop orange", "bop purple", "magma",  "gray"]
    superscript_numbers = [ "\u2070", "\u00b9","\u00b2","\u00b3","\u2074","\u2075","\u2076","\u2077","\u2078", "\u2079" ]
    illegal_columns = ['ID', 'ID03', 'mp', 'µm', 'Time in Seconds '] + ['TB'] # TB contains some image data - but I do not know what exactly

    def get_first_element(self):
        return list(self.elements.keys())[0]

    # find elements and prepare element arrays
    def process_elements(self, line_data):
        if len(self.elements) == 0:
            columns = line_data.columns.tolist()
            for c in self.illegal_columns:
                if c in columns: columns.remove(c)
            for p, element in enumerate(columns):
                e = re.split('(\d+)',element)
                if len(e) == 3:
                    i = ""
                    for n in e[1]:
                        i += self.superscript_numbers[int(n)]
                    self.elements[element] = i+e[0]
                else:
                    self.elements[element] = element

    # interpolate lines in x-direction
    def strech_img(self, img, stretch_x=None):

        if stretch_x is None: stretch_x = self.settings["stretch_x"]

        _, j = img.shape
        line_a = img[:,0]
        new_img = img
        for pos in range(j-1):
            line_b = img[:,pos+1]

            #linear interpolation
            if stretch_x == 1:
                new_img = np.insert(new_img, 2*pos+1, ((line_a + line_b) / 2), axis=1)
            else:
                for i in range(stretch_x):
                    new_img = np.insert(new_img, (stretch_x+1)*pos+(i+1), ((line_b - line_a) / (stretch_x+1) * (i+1) + line_a), axis=1)

            line_a = line_b

        return new_img

    # change data format, the orientation of the image and remove 1% outliers
    def optimize_img(self, img):
        # in np array umwandeln und um 90° drehen
        np_img = np.rot90( np.flip( np.array(img), 0 ), 3 )
        element_max = np.percentile(np_img, 99)
        #plt.hist(np_img.flatten(), bins = range(0,round(element_max[element]), 15))
        #plt.title("histogram {}".format(element))
        #plt.show()
        np_img = np.clip(np_img, 0, element_max)  / element_max
        if self.settings["stretch_x"] > 0: np_img = self.strech_img( np_img )

        return np_img, element_max # normieren

    def process_ls_icp_ms_line(self, line_data):
        self.process_elements( line_data )

        for element in self.elements.keys():
            if not element in self.images: self.images[element] = []
            if settings["smooth_y"] > 0:
                self.images[element].append( gaussian_filter1d(line_data[element].tolist(), self.settings["smooth_y"]) )
            else:
                self.images[element].append( line_data[element].tolist() )

    def load_excel( self ):
        xl = pd.ExcelFile(self.settings["excel_file"])
        sheets = xl.sheet_names

        assert len(sheets) > 1, "The excel file is not formatted as expected! Every data line has to be in a individual excel sheet."

        sheets.sort()
        for sheet in sheets:
            line_data = xl.parse(sheet)
            self.process_ls_icp_ms_line(line_data)

    def load_raw( self ):
        cnt = 0

        for file in os.listdir(self.settings["workingDirectory"]):
            filename = os.fsdecode(file)
            if not '~lock' in filename and '.xl' in filename:
                line_data  = pd.read_csv(self.settings["workingDirectory"] + filename, header=1)
                self.process_ls_icp_ms_line(line_data)
                cnt += 1

        assert cnt != 0, "Did not find any *.xl file in '{}' !".format(self.settings["workingDirectory"])

    def save_image_by_element( self, element, output_directory = None):
        if output_directory is None: output_directory = self.settings["outputDirectory"]
        assert output_directory != '', "no output directory given!"

        img, _ = self.optimize_img( self.images[element] )

        # convert to 16 bit TIF
        img *= 65535
        img = img.astype(np.uint16)

        # save scaling for programs like imagej
        info = {}
        info[282] = round(1/(7/13*1000), 6)
        info[283] = round(1/(self.spot_distance_y*1000), 6)
        info[270] = "ImageJ=" + 'LA_ICP_MS' + "\nunit=" + 'nm'

        pil_img = Image.fromarray(img.astype(np.uint16))
        pil_img.save( output_directory + 'LA-ICP-MS_{}.tif'.format(element), tiffinfo = info )

    def save_images( self ):
        if self.settings["outputDirectory"] != '':
            if not os.path.isdir( self.settings["outputDirectory"] ):
                os.mkdir( self.settings["outputDirectory"] )
            if ( self.verbose ) : print( "Images will be stored in : "  + self.settings["outputDirectory"] )
            for element in self.elements.keys():
                self.save_image_by_element( element )
        else:
            if ( self.verbose ) : print( "Images won't be stored since no output directory is given!" )

    def show_single_image( self, element=None ):
        if element is None: element = self.get_first_element()
        img, _ = self.optimize_img( self.images[element] )
        plt.rcParams['figure.figsize'] = [12, 12]
        plt.title('results for {}'.format(self.elements[element]))
        plt.imshow(img, aspect=self.spot_distance_y/self.spot_distance_x, cmap='gray', interpolation=None)

    def pre_processed_images(self):
        if len(self.np_images) == 0:
            for element in self.elements.keys():
                self.np_images[element], self.element_max[element] = self.optimize_img( self.images[element] )

    # selected_elements has to be a list of elements contained in the data.
    # e.g.: ['Na23', 'Mg24', 'Al27', 'K39']
    def show_image_set( self, selected_elements = [] ):
        np_images = {}
        element_max = {}

        for element in selected_elements:
            if not element in self.elements.keys():
                print('WARNING: {} does not exist in the dataset!')
                selected_elements.remove(element)

        if len(selected_elements) == 0: selected_elements = self.elements.keys()

        self.pre_processed_images()

        with napari.gui_qt():
            viewer = napari.Viewer()
            for i, element in enumerate(selected_elements):
                new_layer = viewer.add_image(self.np_images[element], name='LA-ICP-MS [{}]'.format(self.elements[element]), scale=self.scaling, colormap=self.colormaps[i], opacity=1/len(selected_elements), blending="additive", rendering="iso")
            viewer.scale_bar.visible = True

    def __init__( self, settings ):
        # process input variables
        self.settings = getBaseSettings()
        for key in settings:
            if key in self.settings:
                self.settings[key] = settings[key]
            else:
                print('settings key "{}" is unknown and will be ignored'.format(key))
        self.verbose = settings["showDebuggingOutput"]

        # get / process filepaths
        if self.settings["workingDirectory"] == '':
            #remove root windows
            root = tk.Tk()
            root.withdraw()
            root.attributes("-topmost", True)
            if self.settings["load_raw"]:
                self.settings["workingDirectory"] = filedialog.askdirectory(title='Please select the directory containing the raw data')
            else:
                self.settings["excel_file"] = filedialog.askopenfilename(initialdir=os.getcwd(), title='Please select the excel file containing LA-ICP-MS data', filetypes=[('Excel-file', '*.xlsx')])
                self.settings["workingDirectory"] = os.path.dirname( self.settings["excel_file"] ) + os.sep

            if self.settings["outputDirectory"] != '':
                self.settings["outputDirectory"] = self.settings["workingDirectory"] + self.settings["outputDirectory"] + os.sep

        if ( self.verbose ) : print( "Selected working directory: " + self.settings["workingDirectory"] )

        if os.path.isdir( self.settings["workingDirectory"] ) :
            self.spot_distance_y = self.settings["spot_distance_y"]
            self.spot_distance_x = 7/(self.settings["stretch_x"] + 1)
            # set basic variables containing all elements and their respective colors in the napari editor
            self.scaling = (self.spot_distance_y, self.spot_distance_x)

            # load data
            if self.settings["excel_file"] != '':
                self.load_excel()
            else:
                self.load_raw()

            if self.verbose:
                first_e = self.get_first_element()
                print('loaded a dataset with the dimensions of {} x {} datapoints and {} elements:'.format( len(self.images[first_e]), len(self.images[first_e][0]), len(self.images) ) )
                print(list(self.elements.values()))

        else:
            assert False, 'directory "{}" does not exist'.format( self.settings["workingDirectory"] )



### actual program start
if __name__ == '__main__':

    programInfo()
    settings = processArguments()

    if ( settings["showDebuggingOutput"] ) : print( "I am living in '{}'".format(settings["home_dir"] ) )

    laser_data = LA_ICP_MS_LOADER(settings)
    laser_data.save_images()
    laser_data.show_image_set()

    print( "Script DONE!" )