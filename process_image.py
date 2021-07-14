# -*- coding: utf-8 -*-

from PIL import Image
import os, sys, re, napari, getopt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from skimage.filters import sobel, threshold_multiotsu
from skimage.segmentation import slic, watershed, mark_boundaries
from skimage import color
import tkinter as tk
from tkinter import filedialog
from molmass import Formula as form

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
        # extracted from a SEM image
        "spot_distance_x"      : 7.0,
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

# returns the mass portion of the relevant element
# e.g.: Al2O3 -> 2*m[Al] / m[Al2O3]
def get_oxide_portion( formula ):
    oxide = form(formula)
    for symbol in oxide._elements:
        if symbol != 'O':
            return form(symbol).mass*oxide._elements[symbol][0] / oxide.mass
            break

# extracting only the element symbol
def get_element_from_isotope( isotope ):
    return re.split('(\d+)', isotope)[0]

class LA_ICP_MS_LOADER:
    elements    = {} # list of elements in the raw data
    raw_image   = {} # raw image data
    images      = {} # smoothed signal data
    np_images   = {} # data as np images with a value range from 0-1
    cal_img_ppm = {} # image data calibrated as ppm
    cal_img_mpo = {} # image data calibrated as m.-% oxide
    element_max = {} # max value in the data
    unit = 'µm'      # unit of the pixel dimensions

    # known oxides
    oxide_dict = {
        'Na': 'Na2O',
        'Mg': 'MgO',
        'Al': 'Al2O3',
        'K' : 'K2O',
        'Ti': 'TiO2',
        'V' : 'V2O5',
        'Cr': 'Cr2O3',
        'Mn': 'MnO2',
        'Zn': 'ZnO',
        'Rb': 'Rb2O',
        'Ba': 'BaO',
        'Ca': 'CaO',
        'Sr': 'SrO',
        'P' : 'P4O6',
        'Cu': 'CuO',
        'Ni': 'NiO',
        'Pb': 'PbO',
        'As': 'As2O3'
    }

    illegal_columns = ['ID', 'ID03', 'mp', 'µm', 'Time in Seconds '] + ['TB'] # TB contains some image data - but I do not know what exactly
    colormaps_napari    = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', "bop blue", "bop orange", "bop purple"]
    colormaps           = [(255,0,0) , (0,255,0), (0,0,255), (0,255,255), (255,0,255), (255,255,0), (255, 204, 18), (18, 255, 247), (0, 255, 145)]
    superscript_numbers = [ "\u2070", "\u00b9","\u00b2","\u00b3","\u2074","\u2075","\u2076","\u2077","\u2078", "\u2079" ]
    subscript_numbers   = [ "\u2080", "\u2081","\u2082","\u2083","\u2084","\u2085","\u2086","\u2087","\u2088", "\u2089" ]

    def get_first_element(self):
        return list(self.elements.keys())[0]

    # find elements and prepare element arrays
    def process_elements(self, line_data):
        if len(self.elements) == 0:
            columns = line_data.columns.tolist()
            for c in self.illegal_columns:
                if c in columns: columns.remove(c)
            for p, element in enumerate(columns):
                e = re.split('(\d+)', element) # get_element_from_isotope does not work here.
                if len(e) == 3:
                    i = ""
                    for n in e[1]:
                        i += self.superscript_numbers[int(n)]
                    self.elements[element] = i+e[0]
                else:
                    self.elements[element] = element

    def get_oxide_conc_factor( self, isotope ):
        # Masseanteil am Oxide
        # zb Anteil von Na an Na2O
        selected_oxide = self.oxide_dict[ get_element_from_isotope( isotope ) ]
        return get_oxide_portion(selected_oxide)*self.cal_dict[isotope]*1000000

    # use the calibration matrix to get an 2D array with calibrated values for:
    #  mpo = mass-%-oxide
    #  ppm = parts per million
    def use_calibration( self, value, element, ppm_mpo='ppm' ):
        if self.cal_dict is None: self.set_calibration_dictionary()
        if ppm_mpo == 'ppm':
            calibrated_value = value/self.cal_dict[element]
            #calibrated_value = round( calibrated_value )
        else:
            calibrated_value = (value/self.get_oxide_conc_factor( element ))*100

        return calibrated_value

    # set the calibration dictionary
    def set_calibration_dictionary( self, cal_dict = None ):
        if not isinstance(cal_dict, dict):
            print("\"cal_dict\" is expected to be an dictionary. E.g.: \{'Na': 1.337, 'K': 0.123\}")
        if cal_dict is None or not isinstance(cal_dict, dict):
            print('Warning: Data is not calibrated! Please add a calibration dictionary to process the correct chemical composition.')
            # filling the dictionary with dummy data
            cal_dict = {}
            for isotope in self.elements.keys():
                cal_dict[isotope] = 1.0
        else:
            for isotope in self.elements.keys():
                if not isotope in cal_dict:
                    print('Warning: {} is missing in cal_dict but exists in the raw data!'.format(self.elements[isotope]))
                    cal_dict[isotope] = 1.0

        for isotope in cal_dict:
            if not get_element_from_isotope( isotope ) in self.oxide_dict:
                print( 'Warning: The oxide of {} is not known! Please edit "process_image.py" first and add the oxide to the "oxide_dict" dictionary.' )
                cal_dict.pop( isotope )

        self.cal_dict = cal_dict


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
    def optimize_img(self, img, remove_outliers=True):
        # in np array umwandeln und um 90° drehen
        np_img = np.rot90( np.flip( np.array(img), 0 ), 3 )
        element_max = np.percentile(np_img, 99)
        if remove_outliers:
            np_img = np.clip(np_img, 0, element_max)  / element_max

        #plt.hist(np_img.flatten(), bins = range(0,round(element_max[element]), 15))
        #plt.title("histogram {}".format(element))
        #plt.show()

        if self.settings["stretch_x"] > 0: np_img = self.strech_img( np_img )

        return np_img, element_max # normieren

    def process_ls_icp_ms_line(self, line_data):
        self.process_elements( line_data )

        for element in self.elements.keys():
            if not element in self.images:
                self.images[element]    = []
                self.raw_image[element] = []
            if self.settings["smooth_y"] > 0:
                self.images[element].append( gaussian_filter1d(line_data[element].tolist(), self.settings["smooth_y"]) )
            else:
                self.images[element].append( line_data[element].tolist() )
            self.raw_image[element].append( line_data[element].tolist() )

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
        return img

    def pre_processed_images(self):
        if len(self.np_images) == 0:
            for element in self.elements.keys():
                self.np_images[element], self.element_max[element] = self.optimize_img( self.images[element] )

    def get_color_by_element(self, element, napari_cmap=True):
        i = list(self.elements.keys()).index(element)

        if len(self.colormaps) < len(self.elements):
            for i in range( len(self.elements)-len(self.colormaps) ):
                self.colormaps.append((255,255,255))
                self.colormaps_napari.append('gray')

        if napari_cmap:
            cmap = self.colormaps_napari[i]
        else:
            c = (self.colormaps[i][0]/255, self.colormaps[i][1]/255, self.colormaps[i][2]/255)
            colors = [(0, 0, 0), c] # first color is black, last is color (r, g, b)
            cmap = LinearSegmentedColormap.from_list(element, colors, N=256)

        return cmap

    def get_calibrated_images(self, element, image=None):
        if image is None:
            image = self.np_images[element]
        if len(self.cal_img_ppm) == 0:
            for e in self.elements.keys():
                max_el = self.element_max[element]
                self.cal_img_ppm[e] = image * self.use_calibration( self.element_max[e], e, 'ppm' )
                self.cal_img_mpo[e] = image * self.use_calibration( self.element_max[e], e, 'mpo' )
        return self.cal_img_ppm[element], self.cal_img_mpo[element]

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

        #with napari.gui_qt():
        viewer = napari.Viewer()
        for i, element in enumerate(selected_elements):
            new_layer = viewer.add_image(self.np_images[element], name='LA-ICP-MS [{}]'.format(self.elements[element]), scale=self.scaling, colormap=self.get_color_by_element(element), opacity=1/len(selected_elements), blending="additive", rendering="iso")
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
            self.spot_distance_x = self.settings["spot_distance_x"]/(self.settings["stretch_x"] + 1)
            self.spot_distance_y = self.settings["spot_distance_y"]
            # set basic variables containing all elements and their respective colors in the napari editor
            self.scaling = (self.spot_distance_y, self.spot_distance_x)

            #self.pixel_scaling['x'] = self.spot_distance_x
            #self.pixel_scaling['y'] = self.spot_distance_y
            self.pixel_area = self.spot_distance_y * self.spot_distance_x

            # load data
            if self.settings["excel_file"] != '':
                self.load_excel()
            else:
                self.load_raw()

            first_e = self.get_first_element()

            #print(laser_data.np_images[element].shape[0] * laser_data.spot_distance_y  )
            #print(laser_data.np_images[element].shape[1] * laser_data.spot_distance_x  )

            self.image_dimensions = (len(self.images[first_e][0])*self.spot_distance_y, (len(self.images[first_e])*(self.settings["stretch_x"]+1)-self.settings["stretch_x"])*self.spot_distance_x, self.unit)

            if self.verbose:
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