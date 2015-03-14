import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from astropy.io import fits

import utils.peakdet
import utils.utils
from scipy.signal import savgol_filter
from sklearn.preprocessing import scale

import sys

def get_data_from_fits(file_path):
    hdu_list = fits.open(file_path)
    data = np.array(hdu_list[0].data)
    hdu_list.close()
    return data

if __name__ != "detect":

    def main(cube_params, file_path, file_path_without_noise, words):


        def detect_lines(cube_params, file_path):

            # Reading cube data
            data = get_data_from_fits(file_path)
            cube_rows = data.shape[1]
            cube_cols = data.shape[2]

            # Extracting cube params
            s_f = cube_params['s_f']
            freq = cube_params['freq']
            spe_bw = cube_params['spe_bw']
            spe_res = cube_params['spe_res']

            # Calculating some other params
            sigma = utils.utils.fwhm2sigma(freq, s_f)
            win_len = int(sigma)
            # The width of the windows must be odd
            if (win_len % 2) == 0:
                win_len -= 1

            # Reading denoise cube for testing purposes
            data_without_noise = get_data_from_fits(file_path_without_noise)

            observed_lines = pd.DataFrame([])

            for xi in xrange(0, cube_rows):
                for yi in xrange(0, cube_cols):

                    if xi == 0 and yi == 0:
                        continue

                    observed_mol = pd.Series(np.zeros([spe_bw]))

                    # The pixel (0, 0) always will be a empty pixel with noise
                    values_noise = data[:,0,0]

                    values_without_noise = data_without_noise[:,xi,yi]
                    values = data[:,xi,yi]

                    # Parameter to reduce noise
                    poli_order = 3
                    mean_noise = np.mean(savgol_filter(values_noise,
                                                        win_len, poli_order))
                    std_noise = np.std(savgol_filter(values_noise,
                                                        win_len, poli_order))

                    # Apllying a Savitzky-Golay first derivative
                    # seven-point averaging algorithm is applied to the raw dataset
                    values_denoised = savgol_filter(values, win_len, poli_order)

                    # Delta of temperature to detect a peak candidate
                    scan_sensibity = 3*std_noise
                    # Array with max and min detected in the spectra
                    maxtab, mintab = utils.peakdet.peakdet(values_denoised, scan_sensibity)

                    # Plot detected max
                    if xi == 1 and yi == 1:
                        plt.title(str(xi) + " " + str(yi))
                        plt.plot(values, color='g', label='Observed')
                        plt.plot(values_without_noise, color='y', label='Without Noise')
                        plt.plot(values_denoised, color='r', label='Observed Denoised')
                        plt.legend(loc='upper right')


                    # while len(maxtab) > 0:
                    for max_line_temp in maxtab[:,1]:
                        # max_line_temp = max(maxtab[:,1])

                        if max_line_temp > mean_noise + scan_sensibity:

                            max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]

                            # Plot max point
                            if xi == 1 and yi == 1:
                                plt.plot(max_line_freq, max_line_temp, 'bs', label='Maxima')

                            # gaussian_fitted = (max_line_temp)*\
                            #                utils.utils.gaussian(
                            #                    np.arange(0, spe_bw, spe_res),
                            #                                     max_line_freq, s_f)

                            # for i in xrange(0, len(values_denoised)):  # TO DO: Bad python
                            #     values_denoised[i] = max(values_denoised[i] - gaussian_fitted[i], mean_noise)

                            observed_mol[int(max_line_freq)] = 1

                            # maxtab, mintab = utils.peakdet.peakdet(values_denoised,scan_sensibity)

                        # else:
                        #     break

                    if xi == 1 and yi == 1:
                        plt.show()

                    observed_mol.index = np.arange(freq - int(spe_bw/2),freq + int(spe_bw/2), spe_res)
                    observed_lines[xi, yi] = observed_mol

            return observed_lines

        def near_obs_prob(freq_theo, X_detected, s_f):
            # Pixel of the cube to train the algorithm
            default_pixel_trainer = (1, 1)

            max_prob = 0
            for freq_obs in X_detected.index:

                if X_detected[default_pixel_trainer].loc[freq_obs] != 0:

                    prob = utils.utils.gaussian(freq_theo, freq_obs, (s_f))
                    if prob == 1:
                        return [prob, freq_obs]
                    elif max_prob < prob:
                        max_prob = prob
                        freq_max_prob = freq_obs
            return [max_prob, freq_max_prob]

        def recal_words(words, X_detected, s_f):

            words_recal = pd.DataFrame(np.zeros(words.shape))
            words_recal.index = words.index
            words_recal.columns = words.columns

            for mol in words.columns:
                # Auxiliar array that holds modified observed frequencies
                changed_prob_freqs = np.array([])

                # The word must have theoretical lines in the window
                if np.mean(words[mol]) == 0:
                    continue

                # Then, those lines in certain frequencies will be reeplaced
                # by the probability of the theoretical line to the nearest
                # observed line (Gaussian decay distance)
                for freq_theo in words[mol].index:
                    if words[mol].loc[freq_theo] != 0:

                        max_prob, freq_obs = near_obs_prob(freq_theo, X_detected,
                                                                s_f)
                        # In order to reeplace the highest probability for each
                        # observed frecuency, we save the freq that had been
                        # reeplaced already.

                        if freq_obs not in changed_prob_freqs:
                            words_recal[mol].loc[freq_obs] = max_prob
                            np.append(changed_prob_freqs, freq_obs)
                        elif max_prob > words_recal[mol].loc[freq_obs]:
                            words_recal[mol].loc[freq_obs] = max_prob
            return words_recal


        X_detected = detect_lines(cube_params, file_path)
        words_recal = recal_words(words, X_detected, cube_params['s_f'])
        # return words, X_detected
        return words_recal, X_detected
