import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from astropy.io import fits

import utils.peakdet
import utils.utils
from scipy.signal import savgol_filter

import sys


if __name__ != "detect":

    def main(cube_params, file_path, words):


        def detect_lines(cube_params, file_path):

            hdu_list = fits.open(file_path)
            data = hdu_list[0].data
            cube_rows = data.shape[1]
            cube_cols = data.shape[2]

            s_f = cube_params['s_f']
            freq = cube_params['freq']
            spe_bw = cube_params['spe_bw']
            spe_res = cube_params['spe_res']

            sigma = utils.utils.fwhm2sigma(freq, s_f)
            win_len = int(sigma)
            # The width of the windows must be odd
            if (win_len % 2) == 0:
                win_len -= 1

            observed_lines = pd.DataFrame([])



            for xi in xrange(1, cube_rows):
                for yi in xrange(1, cube_cols):

                    observed_mol = pd.Series(np.zeros([spe_bw]))

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
                            gaussian_fitted = (max_line_temp)*\
                                              utils.utils.gaussian(
                                                  np.arange(0, spe_bw, spe_res),
                                                                   max_line_freq, s_f)

                            for i in xrange(0, len(values_denoised)):  # TO DO: Bad python
                                values_denoised[i] = max(values_denoised[i] - gaussian_fitted[i], mean_noise)

                            observed_mol[int(max_line_freq), ] = 1

                            maxtab, mintab = utils.peakdet.peakdet(values_denoised,scan_sensibity)

                        else:
                            break

                    observed_mol.index = np.arange(freq - int(spe_bw/2),freq + int(spe_bw/2), spe_res)
                    observed_lines[xi, yi] = observed_mol

            return observed_lines

        def near_word_distance(freq_detected, word, s_f):
            if np.mean(word) == 0:
              return 0

            max_distance = 0
            for freq_word in word.index:

                if word[freq_word] != 0:

                    distance = utils.utils.gaussian(freq_word, freq_detected, s_f)

                    if max_distance < distance:
                        max_distance = distance

            return max_distance

        def recalculate_words(words, X_detected, s_f):

            # Pixel of the cube to train the algorithm
            default_pixel_trainer = (1L, 1L)

            words_recalculated = pd.DataFrame(np.zeros(words.shape))
            words_recalculated.index = words.index
            words_recalculated.columns = words.columns

            for freq_detected in X_detected.index:
                # If the observed spectra has a observed line
                # The theoretical words must have lines in the window
                # Then, those lines in certain frequencies will be reeplaced
                # by the distance of the theoretical line to the nearest
                # observed line (Gaussian decay distance)
                if X_detected[default_pixel_trainer].loc[freq_detected] != 0:

                      for mol in words.columns:
                          words_recalculated[mol].loc[freq_detected] = near_word_distance \
                                          (
                                            freq_detected,
                                            words[mol],
                                            s_f
                                          )
            return words_recalculated


        X_detected = detect_lines(cube_params, file_path)
        words_recalculated = recalculate_words(words, X_detected, cube_params['s_f'])
        # return words, X_detected
        return words_recalculated, X_detected


