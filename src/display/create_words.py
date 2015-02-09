from asidopy import *
import math
import sys
import pandas as pd

# Without redshift (Rvel = 0)
# Temp 300 Kelvin
rvel = 0.0
temp = 300.0

if __name__ != "create_words":

    # Function to create the words necessary to fit a sparse coding model
    # to the observed spectra in the previous created cube. It uses:
    #
    #         - freq    : spectral center (frequency)
    #         - spe_res : spectral resolution
    #         - spe_bw  : spectral bandwidth
    #         - s_f     : the width of the spectral lines (fwhm)
    # Returns a DataFrame with a vector for each theoretical line for each isotope
    # in molist
    def gen_words(isolist, cube_params, cube_name):

        log = sys.stdout
        dbpath = 'ASYDO'

        dictionary = pd.DataFrame([])

        for iso in cube_params['molist']:

            univ=vu.Universe(log)
            univ.create_source('word-'+ iso, 0.0, 0.0)
            s_x=1
            s_y=1
            rot=0
            s_f=cube_params['s_f']
            angle=math.pi
            model=vu.IMCM(log,dbpath,iso,temp,('normal',s_x,s_y,angle),
                                              ('skew',s_f,0),
                                              ('linear',angle,rot))
            model.set_radial_velocity(rvel)
            univ.add_component('word-'+ iso, model)

            alpha = 0.0
            delta = 0.0

            line = univ.gen_cube(
                                'observerd', alpha, delta, freq, 10, 20,
                                 cube_params['spe_res'], cube_params['spe_bw']
                                )

            freq = [int(double) for double in cube.freq_axis]
            values = cube.get_spectrum()

            dictionary[iso] = values
            dictionary.index = freq

            plt.plot(cube.freq_axis, cube.get_spectrum())
            plt.show()

        # range = str(int((freq - spe_bw/2.0)/1000.0)) + " - " + str(int((freq + spe_bw/2.0)/1000.0))
        # dictionary.T.to_csv("dictionary/" + range + ".csv")
        return dictionary

if __name__ == "__main__":
    pass
