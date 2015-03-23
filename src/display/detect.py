import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from astropy.io import fits

import utils.peakdet
from scipy.signal import savgol_filter

import sys

# ## Helper constants ###
SPEED_OF_LIGHT = 299792458.0
S_FACTOR = 2.354820045031  # sqrt(8*ln2)
KILO = 1000
DEG2ARCSEC = 3600.0

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

    def __init__(self, cube_params, file_path):

        # Parameter to reduce noise
        self.poli_order = 3
        # Cube path
        self.file_path = file_path
        # Cube parameters
        self.freq = cube_params["freq"]
        self.spe_bw = cube_params["spe_bw"]
        self.spe_res = cube_params["spe_res"]
        self.s_f = cube_params["s_f"]
        self.mean_noise, self.std_noise = self.get_noise_parameters_from_fits()
        self.threshold = self.get_thresold_parameter()

    def get_data_from_fits(self):
        hdu_list = fits.open(self.file_path)
        data = np.array(hdu_list[0].data)
        hdu_list.close()
        return data

    def get_noise_parameters_from_fits(self):
        win_len = self.get_win_len_from_s_f()

        # The pixel (0, 0) always will be a empty pixel with noise
        values_noise = self.get_data_from_fits()[:,0,0]
        mean_noise = np.mean(savgol_filter(values_noise,
                                            win_len,
                                            self.poli_order))
        std_noise = np.std(savgol_filter(values_noise,
                                            win_len,
                                            self.poli_order))
        return mean_noise, std_noise


    def get_thresold_parameter(self):
        return self.mean_noise + 3*self.std_noise


    def get_win_len_from_s_f(self):

        # Calculating some other params
        sigma = fwhm2sigma(self.freq, self.s_f)
        win_len = int(sigma)
        # The width of the windows must be odd
        if (win_len % 2) == 0:
            win_len -= 1
        return win_len

    def get_freq_index_from_params(self):
        return np.arange(self.freq - int(self.spe_bw/2),
                         self.freq + int(self.spe_bw/2),
                         self.spe_res)

    def detect_lines_simple(self, maxtab, values_denoised):

        observed_lines = pd.Series(np.zeros([self.spe_bw]))
        for max_line_temp in maxtab[:,1]:

            if max_line_temp > self.threshold:

                max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]

                # Plot max point
                plt.plot(max_line_freq, max_line_temp, 'bs')

                observed_lines.iloc[int(max_line_freq)] = 1
        return observed_lines

    def detect_lines_subtracting_gaussians(self, maxtab, values_denoised):

        observed_lines = pd.Series(np.zeros([self.spe_bw]))
        while len(maxtab) > 0:
            max_line_temp = max(maxtab[:,1])

            if max_line_temp > self.threshold:

                max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]

                # Plot max point
                plt.plot(max_line_freq, max_line_temp, 'bs')

                gauss_fitt = (max_line_temp)*\
                            gaussian(np.arange(0, self.spe_bw,  self.spe_res),
                                     max_line_freq,  self.s_f)


                for i in xrange(0, len(values_denoised)):
                    values_denoised[i] = max(values_denoised[i] - gauss_fitt[i],
                                             self.mean_noise)

                observed_lines.iloc[int(max_line_freq)] = 1
                maxtab, mintab = utils.peakdet.peakdet(values_denoised,
                                                       self.threshold)
            else:
                break
        return observed_lines

    def get_lines_from_fits(self):
        lines = pd.Series([dict() for i in range(self.spe_bw)])
        lines.index = self.get_freq_index_from_params()

        i = 3
        hdu_list = fits.open(self.file_path)
        while(i < len(hdu_list)):
            for line in hdu_list[i].data:
                """
                    line[1] : Formula
                    line[3] : Frequency (MHz)
                    line[6] : Temperature (No unit)
                """
                lines[int(line[3])][line[1]] = line[6]
            i = i + 3
        hdu_list.close()
        return lines

    def detect_lines(self, train_pixel):
        # Reading
        #
        # Reading cube data
        data = self.get_data_from_fits()

        # Reading denoise cube for testing purposes
        # data_without_noise = get_data_without_noise_from_fits()

        # Pre-processing
        #
        # The training pixel values
        # values_without_noise = data_without_noise[:,xi,yi]
        values = data[:, train_pixel[0], train_pixel[1]]

        # Apllying a Savitzky-Golay first derivative
        # seven-point averaging algorithm is applied to the raw dataset
        win_len = self.get_win_len_from_s_f()
        values_denoised = savgol_filter(values, win_len, self.poli_order)

        # Detecting
        #
        # Plot detected max
        plt.plot(values, color='g', label='Observed')
        plt.xlabel('Relative Frequency [MHz]')
        plt.ylabel('Temperature [No unit]')
        # plt.plot(values_without_noise, color='y', label='Without Noise')
        plt.plot(values_denoised, color='r')#, label='Observed Denoised')
        # plt.legend(loc='upper right')
        # Array with max and min detected in the spectra
        maxtab, mintab = utils.peakdet.peakdet(values_denoised, self.threshold)
        observed_lines = self.detect_lines_subtracting_gaussians(maxtab, values_denoised)
        observed_lines.index = self.get_freq_index_from_params()
        plt.axhspan(0, self.threshold, facecolor='0.5')
        plt.show()

        return observed_lines

    def near_obs_prob(self, freq_theo, X_detected, s_f):
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

    def recal_words(self, words, X_detected):

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

                    max_prob, freq_obs = self.near_obs_prob(freq_theo,
                                                       X_detected,
                                                       self.s_f)
                    # In order to reeplace the highest probability for each
                    # observed frecuency, we save the freq that had been
                    # reeplaced already.
                    if freq_obs not in changed_prob_freqs:
                        words_recal[mol].loc[freq_obs] = max_prob
                        np.append(changed_prob_freqs, freq_obs)
                    elif max_prob > words_recal[mol].loc[freq_obs]:
                        words_recal[mol].loc[freq_obs] = max_prob
        return words_recal

    def train(self, train_pixel, words):
        X_detected = self.detect_lines(train_pixel)
        words_recal = self.recal_words(words, X_detected)
        return words_recal, X_detected

    def test(self, test_pixels, words):
        X_detected = detect_lines(test_pixel)
        words_recal = recal_words(words, X_detected)
        return confusion_matrix
