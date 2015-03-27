import numpy as np
import pandas as pd
import math
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

def gaussian_weighted(x, mu, sig, w):
    """

      """
    return np.power(w, 2.) * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def fwhm2sigma(freq,fwhm):
    """
      Compute the sigma in MHz given a frequency in MHz and a fwhm in km/s
      """
    sigma = (fwhm * 1000 / S_FACTOR) * (freq / SPEED_OF_LIGHT)
    return sigma


class Detect:

    def __init__(self, cube_params, file_path, train_pixel):

        # Cube parameters
        self.freq = cube_params["freq"]
        self.spe_bw = cube_params["spe_bw"]
        self.spe_res = cube_params["spe_res"]
        self.s_f = cube_params["s_f"]
        # Parameter to reduce noise
        self.poli_order = 2
        # Cube data
        self.file_path = file_path
        self.data = self.get_data_from_fits()
        # Cube values
        self.values = self.get_values_filtered_normalized(train_pixel)
        self.mean_noise, self.std_noise = self.get_noise_parameters_from_fits()
        self.threshold = self.get_thresold_parameter()

        self.detected_lines = np.zeros([self.spe_bw])
        self.detected_temps = np.zeros([self.spe_bw])
        self.detected_gauss = np.zeros([self.spe_bw])

        self.max_line_freq = []
        self.max_line_temp = []

    def get_data_from_fits(self):
        hdu_list = fits.open(self.file_path)
        data = np.array(hdu_list[0].data)
        hdu_list.close()
        return data

    def get_values_filtered_normalized(self, train_pixel):
        # Pre-processing
        #
        # The training pixel values
        values = self.data[:, train_pixel[0], train_pixel[1]]

        # Apllying a Savitzky-Golay first derivative
        # n#-sigma-point averaging algorithm is applied to the raw dataset
        win_len = self.get_win_len_from_s_f()
        values = savgol_filter(values, win_len, self.poli_order)
        # Normalize by the maximum of the serie
        values = values/np.max(values)
        return values

    def get_noise_parameters_from_fits(self):
        win_len = self.get_win_len_from_s_f()

        # The pixel (0, 0) always will be a empty pixel with noise
        values_noise = self.data[:,0,0]
        values = self.data[:,0,1]
        values = savgol_filter(values, win_len, self.poli_order)

        values_noise = savgol_filter(values_noise, win_len, self.poli_order)
        values_noise = values_noise/np.max(values)
        mean_noise = np.mean(values_noise)
        std_noise = np.std(values_noise)
        return mean_noise, std_noise


    def get_thresold_parameter(self):
        return max(self.mean_noise, 0) + 3*self.std_noise


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

    def plot_data(self):
        values = self.data[:,0,1]
        plt.plot(values, color='r', label='Observed spectra')
        plt.xlabel('Relative Frequency [MHz]')
        plt.ylabel('Temperature [Normalized]')
        plt.legend(loc='upper right')
        plt.show()

    def plot_detected_lines(self):
        # Plot detected max
        plt.plot(self.values, color='r', label='Observed Filtered')
        plt.xlabel('Relative Frequency [MHz]')
        plt.ylabel('Temperature [Normalized]')
        plt.legend(loc='upper right')
        # Plot Array with max and min detected in the spectra
        # Plot max point
        for i in range(len(self.max_line_freq)):
            plt.plot(self.max_line_freq[i], self.max_line_temp[i], 'bs')
        plt.axhspan(0, self.threshold, alpha=0.5)
        plt.show()

    def detect_lines_simple(self, maxtab):

        for max_line_temp in maxtab[:,1]:

            if max_line_temp > max(self.mean_noise, 0):

                max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]

                self.max_line_freq.append(max_line_freq)
                self.max_line_temp.append(max_line_temp)

                # Set 1 as value of the line
                self.detected_lines[int(max_line_freq)] = 1

                # Save the temp for the property
                self.detected_temps[int(max_line_freq)] += max_line_temp

    def detect_lines_subtracting_gaussians(self, maxtab, values):

        val = np.array(values)

        while len(maxtab) > 0:
            max_line_temp = max(maxtab[:,1])

            if max_line_temp > max(self.mean_noise, 0):

                max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]

                self.max_line_freq.append(max_line_freq)
                self.max_line_temp.append(max_line_temp)

                # Fit the gaussian
                gauss_fitt = (max_line_temp)*\
                            gaussian(np.arange(0, self.spe_bw,  self.spe_res),
                                     max_line_freq,  self.s_f)

                # Subtract the gaussian
                for i in range(0, len(values)):
                    val[i] = max(val[i] - gauss_fitt[i], self.mean_noise)

                # Save the observation with gaussians
                self.detected_gauss += gauss_fitt

                # Set 1 as value of the line
                self.detected_lines[int(max_line_freq)] = 1
                # Save the temp for the property
                self.detected_temps[int(max_line_freq)] += max_line_temp

                maxtab, mintab = utils.peakdet.peakdet(val,  self.threshold)
            else:
                break

    def get_lines_from_fits(self):
        #lines = pd.Series([dict() for i in range(self.spe_bw)])
        #lines.index = self.get_freq_index_from_params()
        lines = pd.Series([])

        i = 3
        hdu_list = fits.open(self.file_path)
        while(i < len(hdu_list)):
            for line in hdu_list[i].data:
                """
                    line[1] : Formula
                    line[3] : Frequency (MHz)
                    line[6] : Temperature (No unit)
                """
                lines[int(line[3])]= line[1]
            i = i + 3
        hdu_list.close()
        return lines

    def get_temps_from_fits(self):
        #lines = pd.Series([dict() for i in range(self.spe_bw)])
        #lines.index = self.get_freq_index_from_params()
        lines = pd.DataFrame([])

        i = 3
        hdu_list = fits.open(self.file_path)
        while(i < len(hdu_list)):
            for line in hdu_list[i].data:
                """
                    line[1] : Formula
                    line[3] : Frequency (MHz)
                    line[6] : Temperature (No unit)
                """
                lines.loc[int(line[3])]= [line[1], line[6]]
            i = i + 3
        hdu_list.close()
        return lines

    def detect_lines(self):
        # Detecting
        #
        # Array with max and min detected in the spectra
        maxtab, mintab = utils.peakdet.peakdet(self.values, self.threshold)
        self.detect_lines_subtracting_gaussians(maxtab, self.values)

        return self.detected_gauss

    def near_obs_freq(self, freq_theo):

        min_dist = self.spe_bw
        near_freq = 0

        for freq_obs in range(0, self.spe_bw):

            if self.detected_lines[freq_obs] != 0:

                dist = math.fabs(freq_theo - freq_obs)

                if dist == 0:
                    return freq_obs
                elif min_dist > dist:
                    min_dist = dist
                    near_freq = freq_obs
        return near_freq

    def near_obs_prob(self, freq_theo, near_freq_obs):
        sigma = fwhm2sigma(self.freq, self.s_f)
        gauss_weight = gaussian_weighted(freq_theo, near_freq_obs,
                                      sigma, self.detected_temps[near_freq_obs])
        factor = 2*sigma
        ini = int(round(near_freq_obs - factor))
        end =int(round(near_freq_obs + factor))
        if ini < 0:
            ini = 0
        if end > self.spe_bw:
            end = self.spe_bw
        window = np.arange(ini, end, self.spe_res)
        gauss_fit = gauss_weight*gaussian(window,
                                      near_freq_obs, sigma)
        return gauss_fit, window

    def recal_words(self, words):

        words_recal = pd.DataFrame(np.zeros(words.shape))
        words_recal.index = np.arange(0, self.spe_bw,  self.spe_res)
        words_recal.columns = words.columns

        for mol in words.columns:
            # The theorethical line will be replaced by the max probability of
            # the nearest observed line (Gaussian decay distance weighted)

            for freq_theo in range(0, self.spe_bw):
                if words.iloc[freq_theo][mol] != 0:
                    nof = self.near_obs_freq(freq_theo)
                    gauss_fit, window = self.near_obs_prob(freq_theo, nof)
                    # Reeplace the highest probability for each theoretical line
                    words_recal[mol].iloc[window] = gauss_fit
                    break
        words_recal.index = words.index
        return words_recal

    def train(self, words):
        X_detected = self.detect_lines()
        words_recal = self.recal_words(words)
        return words_recal, X_detected

    def test(self, recal_words):
        X_detected = detect_lines(test_pixel)
        words_recal = recal_words(words)
        return confusion_matrix
