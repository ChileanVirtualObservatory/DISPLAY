import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import display.create_cube
import display.create_words
import display.detect
from sklearn.svm import SVC

def load_dictionary():
    input_file = open('dictionary.pkl', "rb" )
    D = pickle.load( input_file )
    input_file.close()
    return D

def get_fortran_array(input):
    fort_array = np.asfortranarray(np.asmatrix(input))
    fort_array = np.asfortranarray(fort_array, dtype= np.double)
    return fort_array


isolist = set(['HC15Nv=0', 'H13CNv2=1', 'H13CNv=0'])
cube_name = 'observed_cube'
cube_name_without_noise = 'observed_cube_without_noise'
#         Creation of a Data cube that needs the following parameters:
#
#         - freq    : spectral center (frequency)
#         - alpha   : right-ascension center
#         - delta   : declination center
#         - spe_res : spectral resolution
#         - spe_bw  : spectral bandwidth
#         - s_f     : the width of the spectral lines (fwhm)
cube_params = {
'freq'     : 604000,
'alpha'    : 0,
'delta'    : 0,
'spe_bw'   : 4000,
'spe_res'  : 1,
's_f'      : 85
            }
# display.create_cube.gen_cube(isolist, cube_params, cube_name)
# display.create_cube.gen_cube(isolist, cube_params, cube_name_without_noise,
#                           bool_noise= False)

# Function to create the words necessary to fit a sparse coding model
# to the observed spectra in the previous created cube. It uses:
#
#         - freq    : spectral center (frequency)
#         - spe_res : spectral resolution
#         - spe_bw  : spectral bandwidth
#         - s_f     : the width of the spectral lines (fwhm)
# Returns a DataFrame with a vector for each theoretical line for each isotope
# in molist
# dictionary = display.create_words.gen_words(molist, cube_params)
# save_dictionary(D)
dictionary = load_dictionary()

file_path = cube_name + '.fits'
file_path_without_noise = cube_name_without_noise + '.fits'
D, X = display.detect.main(cube_params, file_path, file_path_without_noise, dictionary)

y_train = get_fortran_array(np.asmatrix(X[(1, 1)]))

test_index_cubes = [(i, j) for i in range(2,7) for j in range(2,4)]
y_test = [get_fortran_array(X[(index)]) for index in test_index_cubes]

Dictionary = np.asfortranarray(D, dtype= np.double)

##############################################################################
# Train classifier
#
# For an initial search, a logarithmic grid with basis
# 10 is often helpful. Using a basis of 2, a finer
# tuning can be achieved but at a much higher cost.
C = 10
gamma = 0.1
clf = SVC(C=C, gamma=gamma)
clf.fit(Dictionary.T, dictionary.columns)
