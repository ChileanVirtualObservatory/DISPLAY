from asidopy import *
import math
import sys

import numpy as np
import pandas as pd

# Without redshift (Rvel = 0)
# Temp 300 Kelvin
rvel = 0.0
temp = 300.0

if __name__ != "create_words":

    # Function to create the words necessary to fit a sparse coding model
    # to the observed spectra in the previous created cube. It uses:
    #
    #         - molist  : dictionary of molecules and its respective isotopes
    #         - freq    : spectral center (frequency)
    #         - spe_res : spectral resolution
    #         - spe_bw  : spectral bandwidth
    #         - s_f     : the width of the spectral lines (fwhm)
    # Returns a DataFrame with a vector for each theoretical line for each
    # isotope in molist
    def gen_words(molist, cube_params):

        log = sys.stdout
        dbpath = 'ASYDO'

        s_f = cube_params['s_f']
        freq = cube_params['freq']
        spe_bw = cube_params['spe_bw']
        spe_res = cube_params['spe_res']

        dictionary = pd.DataFrame([])

        for mol in molist:

            for iso in molist[mol]:
                univ=vu.Universe(log)
                univ.create_source('word-'+ iso)
                s_x=1
                s_y=1
                rot=0
                s_f=cube_params['s_f']
                angle=math.pi
                model=vu.IMCM(log,dbpath,iso,temp,('normal',s_x, s_y, angle),
                                                  ('skew', s_f, 0),
                                                  ('linear', angle, rot))
                model.set_radial_velocity(rvel)
                univ.add_component('word-'+ iso, model)

                lines = univ.gen_cube('observerd', freq, spe_res, spe_bw)

                if len(lines.hdulist) > 1:
                    for line in lines.hdulist[1].data:
                        word =  np.array(np.zeros(len(lines.get_spectrum())))
                        '''
                            line[0] : line_code
                            line[1] : relative freq at the window
                        '''
                        word[line[1]] = 1
                        dictionary[line[0]] = word
        dictionary.index = np.arange(freq - int(spe_bw/2),
                                     freq + int(spe_bw/2),
                                     spe_res)

        # range = str(int((freq - spe_bw/2.0)/1000.0)) + " - " + str(int((freq + spe_bw/2.0)/1000.0))
        # dictionary.T.to_csv("dictionary/" + range + ".csv")
        return dictionary

if __name__ == "__main__":
    pass
