from __future__ import division
from sklearn.decomposition import SparseCoder
import numpy as np
from sklearn.metrics import mean_squared_error
import math
import scipy.sparse as ssp
import matplotlib.pyplot as plt
import sys
import spams




if __name__ != "training":

    def normalize(dictionary):
        dictionary = dictionary.T
        for i in xrange(0, dictionary.shape[0]):
            std = np.std(dictionary[i])
            if std != 0:
                dictionary[i] = (dictionary[i] - np.mean(dictionary[i])) / std
        return dictionary.T

    def fit(D, X):

        dictionary = D

        X = np.asfortranarray(np.asmatrix(X).T)
        X = np.asfortranarray(X, dtype= np.float64)

        # Normalization
        # dictionary = normalize(dictionary)

        D = np.asfortranarray(dictionary.T, dtype= np.float64)

        # param = {
        #     'lambda1' : 0.25, # not more than 20 non-zeros coefficients
        #     'L' : 20, # not more than 20 non-zeros coefficients (default: min(size(D,1),size(D,2)))
        #     # parameter of the optimization procedure are chosen
        #     'numThreads' : -1, # number of processors/cores to use; the default choice is -1
        #     # and uses all the cores of the machine
        #     'mode' : 0}        # penalized formulation
        # alpha = spams.lasso(X, D, **param)

        param = {
            'lambda1' : 1000, # not more than 20 non-zeros coefficients
            # 'lambda2' : 0.1,
            # 'L': 1,
            'pos' : True,
            'mode' : 0,
            'ols' : True,
            'numThreads' : -1} # number of processors/cores to use; the default choice is -1
                               # and uses all the cores of the machine
        alpha = spams.lasso(X, D, **param)

        alpha = alpha.toarray()

        for p in xrange(0, len(index)):
          print str(index[p]) + ": " + str(alpha[p])

        # for i in range(0,len(words)):
        #    print str(sys.argv[1] == words[i]) + " " + words[i] + " " + str(alpha[i])

        total = np.inner(D, alpha.T)

        f, axarr = plt.subplots(1, 1)
        # axarr.plot(dictionary[0,:])


        plt.plot(y, color='r', label='Observed')
        plt.plot(total, color='b', label='Recovered')

        plt.legend(loc='upper right')

        plt.show()

        return alpha
