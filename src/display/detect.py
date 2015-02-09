import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits

import utils.peakdet
import utils.utils
from scipy.signal import savgol_filter


if __name__ != "detect":
    # s_f : parameter to detect lines
    def detect_lines(cube_params, file_path):

        hdu_list = fits.open(file_path)
        data = hdu_list[0].data
        cube_rows = data.shape[1]
        cube_cols = data.shape[2]

        # Parameter to reduce noise
        s_f = cube_params['s_f']
        freq = cube_params['freq']
        spe_bw = cube_params['spe_bw']
        spe_res = cube_params['spe_res']

        sigma = utils.utils.fwhm2sigma(freq, s_f)
        win_len = int(sigma)
        # The width of the windows must be odd
        if (win_len % 2) == 0:
            win_len -= 1

        observed_lines = np.zeros([spe_bw, cube_rows, cube_cols])

        for xi in xrange(1, cube_rows):
            for yi in xrange(1, cube_cols):

                # The pixel (0, 0) always will be a empty pixel with noise
                values_noise = data[:,0,0]

                values = data[:,xi,yi]

                # Parameter to reduce noise
                poli_order = 3
                mean_noise = np.mean(savgol_filter(values_noise, win_len, poli_order))
                std_noise = np.std(savgol_filter(values_noise, win_len, poli_order))

                values_denoised = savgol_filter(values, win_len, poli_order)

                # Delta of temperature to detect a peak candidate
                scan_sensibity = 3*std_noise
                # Array with max and min detected in the spectra
                maxtab, mintab = utils.peakdet.peakdet(values_denoised, scan_sensibity)


                while len(maxtab) > 0:
                    max_line_temp = max(maxtab[:,1])

                    if max_line_temp>mean_noise + scan_sensibity:

                        max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]
                        gaussian_fitted = (max_line_temp)*utils.utils.gaussian(np.arange(0,spe_bw,spe_res), max_line_freq, s_f)

                        for i in xrange(0, len(values_denoised)):  # TO DO: Bad python
                            values_denoised[i] = max(values_denoised[i] - gaussian_fitted[i], mean_noise)

                        observed_lines[int(max_line_freq), xi, yi] = 1

                        maxtab, mintab = utils.peakdet.peakdet(values_denoised,scan_sensibity)

                    else:
                        break

                # observed_lines.index = np.arange(freq - int(spe_bw/2),freq + int(spe_bw/2),spe_res)
        return observed_lines

      def recalculate_lines(D, X):

          cube_rows = X.shape[1]
          cube_cols = X.shape[2]
            for xi in xrange(1, cube_rows):
                for yi in xrange(1, cube_cols):
                    for f in xrange(1, len(D)):
                        X[f, xi, yi] = near_word(D, )


