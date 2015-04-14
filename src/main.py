"""
This file is part of ChiVO
Copyright (C) Andres Riveros

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import display.create_cube
import display.create_words

from display.detect import *

import spams

# Dictionary of molecules and its respective isotopes.
molist = {
            'CO' : ('COv=0','COv=1','13COv=0','C18O','C17O','13C17O','13C18O'),
            # Carbon Monoxide

            # 'NH2' : ('NH2'), # Amidogen

            'N2H' : ('N2H+v=0', 'N2D+', '15NNH+', 'N15NH+'), # Diazenylium

            'CN' : ('CNv=0', '13CN', 'C15N'), # Cyanide Radical

            'HCN' : ('HCNv=0', 'HCNv2=1', 'HCNv2=2','HCNv3=1', 'HC15Nv=0',
                     'H13CNv2=1', 'H13CNv=0', 'HCNv1=1', 'DCNv=0',
                     'DCNv2=1', 'HCNv2=4', 'HCNv2=1^1-v2=4^0'),
            # Hydrogen Cyanide

            # 'H2CN' : ('H2CN'), # Methylene amidogen

            'CS' : ('CSv=0', '13C34Sv=0', 'C36Sv=0', 'C34Sv=0', 'CSv=1-0',
                    '13CSv=0', 'C33Sv=0', 'CSv=1', 'C34Sv=1'),
            # Carbon Monosulfide

            'CCS' : ('CCS', 'C13CS', '13CCS', 'CC34S'), # Thioxoethenylidene

            'H2S' : ('H2S', 'H234S', 'D2S'), # Hydrogen sulfide

            'H2CS' : ('H2CS', 'H213CS', 'H2C34S'), # Thioformaldehyde

            'SO2' : ('SO2v=0', '33SO2', '34SO2v=0', 'SO2v2=1', 'OS18O', 'OS17O'),
            # Sulfur Dioxide

            'H2CO' : ('H2CO', 'H2C18O', 'H213CO'), # Formaldehyde

            'HCO' : ('HCO+v=0', 'HC18O+', 'HC17O+', 'H13CO+'), # Formylium

            # 'HC3N' : ('HC3Nv=0'), # Cyanoacetylene

            'HC5N' : ('HC5Nv=0', 'HC5Nv11=1', 'HCC13CCCN', 'HCCCC13CN',
                      'HCCC13CCN', 'H13CCCCCN'), # Cyanobutadiyne

            'CH3OH' : ('CH3OHvt=0', '13CH3OHvt=0 ', 'CH318OH', 'CH3OHvt=1 ',
                       '13CH3OHvt=1 ') # Methanol
          }

def save_dictionary(D):
    output = open('dictionary.pkl', 'wb')
    pickle.dump(D, output)
    output.close()

def load_dictionary():
    input_file = open('dictionary.pkl', "rb" )
    D = pickle.load( input_file )
    input_file.close()
    return D

def get_fortran_array(input):
    fort_array = np.asfortranarray(np.asmatrix(input)).T
    fort_array = np.asfortranarray(fort_array, dtype= np.double)
    return fort_array

def show_alphas_in_isolist(alpha, isolist):
    for p in xrange(0, len(alpha)):
        #if alpha[p] != 0:
        if dictionary.columns[p] in isolist:
            print "*" + dictionary.columns[p] + ": " +  str(alpha[p])
        else:
            print dictionary.columns[p] + ": " +  str(alpha[p])

def show_existent_alphas(alpha):
    for p in xrange(0, len(alpha)):
        if alpha[p] != 0:
            print dictionary.columns[p] + ": " +  str(alpha[p])

def show_words(dictionary, isolist):
    for p in isolist:
        plt.plot(dictionary[p])
    plt.show()

def show_recall(Results):
    for fscore in Results['Recall'].index:
        print fscore + "   " + str(Results['Recall'].loc[fscore])

def show_precision(Results):
    for fscore in Results['Precision'].index:
        print fscore + "   " + str(Results['Precision'].loc[fscore])

def show_fscore(Results):
    for fscore in Results['F-Score'].index:
        print fscore + "   " + str(Results['F-Score'].loc[fscore])

def get_molecule_from_isotope(searched_iso, mydict):
    for mol, iso_list in mydict.items():
        for iso in iso_list:
            if iso == searched_iso:
                return mol

def graph_sparse_coding(Detector, dictionary_recal):
    lines = Detector.get_lines_from_fits()
    # Shows lines really present
    for freq in lines.index:
        plt.axvline(x=freq, ymin=0, ymax= 1, color='g')
        plt.text(freq, 1.0, lines[freq], size='14', rotation='vertical')
    # Show predicted classes
    for freq in range(cube_params['spe_bw']):
        i = 0.0
        if Detector.detected_lines[freq] != 0:
            for mol_ix in range(len(dictionary_recal.columns)):
                line_name = dictionary_recal.columns[mol_ix]
                if dictionary_recal[line_name].iloc[freq] != 0 and alpha[mol_ix] != 0:
                    plt.axvline(x=dictionary_recal.index[freq], ymin=0, ymax= 0, color='g')
                    plt.text(dictionary_recal.index[freq], 0.15 + i, isotope_name, size='14', rotation='vertical')
                    i += 0.3
    xaxis = Detector.get_freq_index_from_params()
    plt.plot(xaxis, y_train, color='r',
             label='Observed')
    plt.plot(xaxis, total, color='b',
             label='Recovered', linestyle='--')
    plt.legend(loc='upper right')
    plt.xlim(xmax = xaxis[-1], xmin = xaxis[0])
    plt.show()


def temporal_test(Detector, dictionary_recal, alpha):
    lines = Detector.get_lines_from_fits()

    set_isotopes = set()
    # Catches and shows lines really present
    for freq in lines.index:
        set_isotopes.add(lines[freq])
    # Catches and shows lines predicted
    for mol_ix in range(0, len(dictionary_recal.columns)):
        if alpha[mol_ix] != 0:
            line_name = dictionary_recal.columns[mol_ix]
            isotope_name = line_name.split('f')[0][:-1]
            set_isotopes.add(isotope_name)

    # Confusion Matrix construction
    MatrixConfusion = pd.DataFrame(np.zeros(
                                            (len(set_isotopes),
                                            len(set_isotopes))
                                            ),
                                    index=set_isotopes,
                                    columns=set_isotopes)

    # print "Observed"
    for freq in range(cube_params['spe_bw']):
        if Detector.detected_lines[freq] != 0:
            # print dictionary_recal.index[freq]
            tot_sum = 0
            TempAlpha = pd.Series([])
            for mol_ix in range(len(dictionary_recal.columns)):
                line_name = dictionary_recal.columns[mol_ix]
                isotope_name = line_name.split('f')[0][:-1]
                if dictionary_recal[line_name].iloc[freq] != 0 and alpha[mol_ix] != 0:
                    TempAlpha[isotope_name] = 1/alpha[mol_ix]
                    tot_sum += 1/alpha[mol_ix]
            TempAlpha = TempAlpha/tot_sum
            # print TempAlpha
            closest_line = lines[min(lines.index, key=lambda x:abs(x-dictionary_recal.index[freq]))]

            for isotope in TempAlpha.index:
                MatrixConfusion[isotope].loc[closest_line] += TempAlpha[isotope]

    Results = pd.DataFrame(np.zeros((len(set_isotopes), 3)),
                      index=set_isotopes, columns=['Precision', 'Recall',
                                                      'F-Score'])

    for isotope in set_isotopes:
        true_positives = 0
        tot = 0
        for row in MatrixConfusion.index:
            #if isotope == row:
            if get_molecule_from_isotope(isotope, molist) == get_molecule_from_isotope(row, molist):
                true_positives += MatrixConfusion.loc[row][isotope]
            tot += MatrixConfusion.loc[row][isotope]
        if tot != 0:
            Results['Precision'].loc[isotope] = true_positives/tot

    for isotope in set_isotopes:
        true_positives = 0
        tot = 0
        for column in MatrixConfusion.columns:
            # if column == isotope:
            if get_molecule_from_isotope(isotope, molist) == get_molecule_from_isotope(column, molist):
                true_positives += MatrixConfusion.loc[isotope][column]
            tot += MatrixConfusion.loc[isotope][column]
        if tot != 0:
            Results['Recall'].loc[isotope] = true_positives/tot

    for isotope in set_isotopes:
        recall = Results['Recall'].loc[isotope]
        precision = Results['Precision'].loc[isotope]
        if recall != 0 or precision != 0:
            Results['F-Score'].loc[isotope] = 2.*(recall*precision)/(recall + precision)
    return Results

# Subset of molist to use in the simulation to generate spectral lines
#
#       @Test: Isotopes with theoretical lines on the Band 9 (602 - 720 Ghz),
#       in a sample of 4 Ghz with resolution 1 Mhz, Sample: [602 - 606]
#       HC15Nv=0
#       H13CNv2=1
#       H13CNv=0
#       H213CS
#       H2C34S
#       SO2v=0
#       33SO2
#       34SO2v=0
#       SO2v2=1
#       OS18O
#       OS17O
#       H2C18O
#       H213CO
#



if __name__ == "__main__":
    # Function to create a fit containing an observed object (a Datacube
    # ALMA-like) using ASYDO Project. Parameters:
    #
    #         - isolist     : list subset of the list of isotopes to generate a cube
    #         - cube_name    : filename of the .fits generated by the simulation
    #         - cube_params : parameters for the simulation of the cube


    isolist = set(['HC15Nv=0', 'H13CNv2=1', 'H13CNv=0', 'H213CS',
                   '34SO2v=0','SO2v2=1','OS18O','OS17O','H2C18O', 'H213CO'])
    cube_name = 'observed_cube'
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
    # display.create_cube.gen_cube(isolist, cube_params, cube_name, white_noise=True)

    # Function to create the words necessary to fit a sparse coding model
    # to the observed spectra in the previous created cube. It uses:
    # Returns a DataFrame with a vector for each theoretical line for each isotope
    # in molist
    # dictionary = display.create_words.gen_words(molist, cube_params)
    # save_dictionary(dictionary)
    dictionary = load_dictionary()

    # Training
    #
    #
    file_path = cube_name + '.fits'
    Detector = Detect(cube_params, file_path, (1, 1))
    dictionary_recal, X = Detector.train(dictionary)

    y_train = get_fortran_array(np.asmatrix(X))
    dictionary_recal_fa = np.asfortranarray(dictionary_recal, dtype= np.double)

    param = {
      'lambda1' : 25, # practically unrestricted = 1000
      # 'lambda2' : 0.1,
      # 'L': 1,
      'pos' : True,
      'mode' : 0,
      'ols' : True,
      'numThreads' : -1} # number of cores to use; the default choice is -1
    # and uses all the cores of the machine
    alpha = spams.lasso(y_train, dictionary_recal_fa, **param).toarray()
    total = np.inner(dictionary_recal_fa, alpha.T)

    # show_alphas_in_isolist(alpha, isolist)

    graph_sparse_coding(Detector, dictionary_recal, alpha)

    # Testing
    #
    # Confusion Matrix construction
    Results = temporal_test(Detector, dictionary_recal, alpha)
