import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from astropy.io import fits

import utils.peakdet
from scipy.signal import savgol_filter

import sys


def gaussian(x, mu, sig):
    """

      """
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def fwhm2sigma(freq,fwhm):
    """
      Compute the sigma in MHz given a frequency in MHz and a fwhm in km/s
      """
    sigma = (fwhm * 1000 / S_FACTOR) * (freq / SPEED_OF_LIGHT)
    return sigma


class Detect:

    def __init__(self, file_path, cube_params):

        # Parameter to reduce noise
        self.poli_order = 3
        self.file_path = file_path
        self.freq = cube_params['freq']
        self.spe_bw = cube_params['spe_bw']
        self.spe_res = cube_params['spe_res']
        self.s_f = cube_params['s_f']

    def get_data_from_fits():
        hdu_list = fits.open(self.file_path)
        data = np.array(hdu_list[0].data)
        hdu_list.close()
        return data

    def get_noise_parameters_from_fits():
        win_len = self.get_win_len_from_s_f()

        # The pixel (0, 0) always will be a empty pixel with noise
        values_noise = self.get_data_from_fits(self.file_path)[:,0,0]
        mean_noise = np.mean(savgol_filter(values_noise,
                                            win_len, poli_order))
        std_noise = np.std(savgol_filter(values_noise,
                                            win_len, poli_order))
        return mean_noise, std_noise

    def get_thresold_parameter(file_path):
        mean_noise, std_noise = self.get_noise_parameters_from_fits()
        return mean_noise + 3*std_noise

    def get_win_len_from_s_f(cube_params):

        # Calculating some other params
        sigma = fwhm2sigma(self.freq, self.s_f)
        win_len = int(sigma)
        # The width of the windows must be odd
        if (win_len % 2) == 0:
            win_len -= 1
        return win_len

    def get_freq_index_from_params():
        return np.arange(freq - int(self.spe_bw/2),
                                         self.freq + int(self.spe_bw/2),
                                         self.spe_res)

    def get_lines_from_fits():
        lines = pd.Series(np.zeros([self.spe_bw]))
        lines.index = get_freq_index_from_params()

        threshold = get_thresold_parameter(self.file_path)
        i = 3
        hdu_list = fits.open(self.file_path)
        while(i < len(hdu_list)):
            for line in hdu_list[i].data:
                if line[6] > threshold:
                    """
                        line[1] : Formula
                        line[3] : Frequency (MHz)
                        line[6] : Temperature (No unit)
                    """
                    lines[line[1]].iloc[int(line[3])] = 1
                i = i + 3
        hdu_list.close()
        return lines

    def detect_lines(default_pixel):
        # Reading
        #
        # Reading cube data
        data = get_data_from_fits(self.file_path)

        # Reading denoise cube for testing purposes
        # data_without_noise = get_data_from_fits(file_path_without_noise)

        # Pre-processing
        #
        # The training pixel values
        # values_without_noise = data_without_noise[:,xi,yi]
        values = data[:,default_pixel['x'], default_pixel['y']]

        # Apllying a Savitzky-Golay first derivative
        # seven-point averaging algorithm is applied to the raw dataset
        win_len = self.get_win_len_from_s_f(self.cube_params)
        values_denoised = savgol_filter(values, win_len, poli_order)

        # Detecting
        #
        # Delta of temperature to detect a peak candidate
        scan_sensibity = get_thresold_parameter(file_path)
        # Array with max and min detected in the spectra
        maxtab, mintab = utils.peakdet.peakdet(values_denoised, scan_sensibity)

        # Plot detected max
        # plt.plot(values, color='b', label='Observed')
        # plt.xlabel('Relative Frequency [MHz]')
        # plt.ylabel('Temperature [No unit]')
        # plt.plot(values_without_noise, color='y', label='Without Noise')
        # plt.plot(values_denoised, color='r')#, label='Observed Denoised')
        # plt.legend(loc='upper right')

        observed_lines = pd.Series(np.zeros([self.spe_bw]))
        # while len(maxtab) > 0:
        for max_line_temp in maxtab[:,1]:
            # max_line_temp = max(maxtab[:,1])

            if max_line_temp > mean_noise + scan_sensibity:

                max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]

                # Plot max point
                # plt.plot(max_line_freq, max_line_temp, 'bs', label='Maxima')

                # gaussian_fitted = (max_line_temp)*\
                #                gaussian(
                #                    np.arange(0, spe_bw, spe_res),
                #                                     max_line_freq, s_f)

                # for i in xrange(0, len(values_denoised)):  # TO DO: Bad python
                #     values_denoised[i] = max(values_denoised[i] - gaussian_fitted[i], mean_noise)

                observed_lines.iloc[int(max_line_freq)] = 1
                # maxtab, mintab = utils.peakdet.peakdet(values_denoised,scan_sensibity)
            # else:
            #     break

        observed_lines.index = get_freq_index_from_params(cube_params)
        return observed_lines

    def near_obs_prob(freq_theo, X_detected, s_f):
        max_prob = 0
        for freq_obs in X_detected.index:

            if X_detected.loc[freq_obs] != 0:

                prob = gaussian(freq_theo, freq_obs, (s_f))
                if prob == 1:
                    return [prob, freq_obs]
                elif max_prob < prob:
                    max_prob = prob
                    freq_max_prob = freq_obs

        return [max_prob, freq_max_prob]

    def recal_words(words, X_detected):

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

                    max_prob, freq_obs = near_obs_prob(freq_theo,
                                                       X_detected,
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

    def train(train_pixel, words):
        X_detected = detect_lines(train_pixel)
        words_recal = recal_words(words, X_detected)
        return words_recal, X_detected

    def test(test_pixels, words):
        X_detected = detect_lines(test_pixel)
        words_recal = recal_words(words, X_detected)
        return confusion_matrix
