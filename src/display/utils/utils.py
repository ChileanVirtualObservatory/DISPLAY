import sys
import math

import numpy as np

SPEED_OF_LIGHT = 299792458.0

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def fwhm2sigma(freq,fwhm):
  """
      Compute the sigma in MHz given a frequency in MHz and a fwhm in km/s
      """
  sigma = (fwhm * 1000 / S_FACTOR) * (freq / SPEED_OF_LIGHT)
  return sigma

def nearest_candidate(potential_line, candidate_list, threshold):
    candidate_distance = sys.maxint
    best_candidate = None
    for candidate in candidate_list:

        distance = math.fabs(candidate.f_max - potential_line)
        if distance < candidate_distance:
            candidate_distance = distance
            best_candidate = candidate

    if not best_candidate:
        return 0

    if candidate_distance < best_candidate.sigma:
        return 1
    return 0

def nearest_candidate_weighted(potential_line, candidate_list, threshold):
    candidate_distance = sys.maxint
    best_candidate = None
    for candidate in candidate_list:

        distance = math.fabs(candidate.f_max - potential_line)
        if distance < candidate_distance:
            candidate_distance = distance
            best_candidate = candidate

    if not best_candidate:
        return 0

    if candidate_distance < best_candidate.sigma:
        # return 2/(best_candidate.sigma*math.sqrt(2*math.pi))*math.exp(-(candidate_distance)**2/(2*best_candidate.sigma**2))
        return np.exp(-np.power(candidate_distance, 2.) / (2 * np.power(best_candidate.sigma, 2.)))
    return 0

def sigma_threshold(path, freq_ini, freq_end):
    sum = 0
    for i in range(freq_ini, freq_end, 2):
        file = 'combined-' + str(i + 1) + '.fits'
        fits = pyfits.open(path + file)
        sum = sum + np.std(fits[0].data[:, 3, 3])

    varsigma = sum/((freq_end-freq_ini)/2-1)
    return varsigma

def mean_threshold(path, freq_ini, freq_end):
    sum = 0
    for i in range(freq_ini, freq_end, 2):
        file = 'combined-' + str(i + 1) + '.fits'
        fits = pyfits.open(path + file)
        sum = sum + np.mean(fits[0].data[:, 3, 3])

    mean = sum/((freq_end-freq_ini)/2-1)
    return mean

def get_cube_from_paths_file(path):
    f = open(path, 'r')

    return [l[:-1] for l in f]

def undo_red_shift(freq, rad_vel):
    return freq*(1.0/(1 + rad_vel*1000.0/SPEED_OF_LIGHT))
