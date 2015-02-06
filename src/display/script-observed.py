import asydopy.vu
import matplotlib.pyplot as plt
import random
import math
import sys
from itertools import combinations
import numpy as np
from scipy.signal import savgol_filter
import utils.peakdet
import utils.utils

# Without redshift (Rvel = 0)
# Temp 300 Kelvin
rvel = 0.0
temp = 300.0

# Molecules and its respective isotopes.
molist = {
            'CO' : ('COv=0','COv=1','13COv=0','C18O','C17O','13C17O','13C18O'),
                    # Carbon Monoxide

            # 'NH2' : ('NH2'), # Amidogen

            'N2H' : ('N2H+v=0', 'N2D+', '15NNH+', 'N15NH+'), # Diazenylium

            'CN' : ('CNv=0', '13CN', 'C15N'), # Cyanide Radical

            'HCN' : ('HCNv=0', 'HCNv2=1', 'HCNv2=2','HCNv3=1', 'HC15Nv=0',
                     'H13CNv2=1', 'H13CNv=0', 'HCNv1=1', 'HCNv3=1', 'DCNv=0',
                     'DCNv2=1', 'HCNv2=4', 'HCNv2=1^1-v2=4^0'),
            # Hydrogen Cyanide

            # 'H2CN' : ('H2CN'), # Methylene amidogen

            'CS' : ('CSv=0', '13C34Sv=0', 'C36Sv=0', 'C34Sv=0', 'CSv=1-0',
                    '13CSv=0', 'C33Sv=0', 'CSv=1', 'C34Sv=1'),
            # Carbon Monosulfide

            'CCS' : ('CCS', 'C13CS', '13CCS', 'CC34S'), # Thioxoethenylidene

            'H2S' : ('H2S', 'H2S', 'H234S', 'D2S'), # Hydrogen sulfide

            'H2CS' : ('H2CS', 'H213CS', 'H2C34S'), # Thioformaldehyde

            'SO2' : ('SO2v=0', '33SO2', '34SO2v=0', 'SO2v2=1'),
            # Sulfur Dioxide

            'OSO' : ('OS18O', 'OS17O'),# Sulfur Dioxide

            'H2CO' : ('H2CO', 'H2C18O', 'H213CO'), # Formaldehyde

            'HCO' : ('HCO+v=0', 'HC18O+', 'HC17O+', 'H13CO+'), # Formylium

            # 'HC3N' : ('HC3Nv=0'), # Cyanoacetylene

            'HC5N' : ('HC5Nv=0', 'HC5Nv11=1', 'HCC13CCCN', 'HCCCC13CN',
                      'HCCC13CCN', 'H13CCCCCN', 'HC13CCCCN'), # Cyanobutadiyne

            'CH3OH' : ('CH3OHvt=0', '13CH3OHvt=0', 'CH318OH', 'CH3OHvt=1',
                       '13CH3OHvt=1') # Methanol
            }

log = sys.stdout
dbpath = 'ASYDO'

"""
list_of_all = set(['H2CS', 'H2CO', 'CH3OH']
isolist = ['COv=0']
Isolist: set of Molecules
for mol in list_of_all:
    for iso in molist[mol]:
        isolist.append(iso)
"""
isolist = set(['H2CS', 'H2CO', 'CH3OH', 'SO2', 'OSO', 'HCN'])

AllPossComb =  sum(map(lambda r: list(combinations(isolist, r)), range(1, len(isolist)+1)), [])

for freq_iteration in range(0, 1):

    for Comb in AllPossComb:

        univ=asydopy.vu.Universe(log)

        if len(Comb) == 1:
            continue

        print Comb

        for key in Comb:

            for mol in molist[key]:

                univ.create_source('observed-'+mol, 0.0, 0.0)
                s_x=random.uniform(50, 150)
                s_y=random.uniform(40, 100)
                rot=random.uniform(10, 150)
                s_f=random.uniform(50, 120)
                angle=random.uniform(0,math.pi)
                model=asydopy.vu.IMCM(log,dbpath,mol,temp,('normal',s_x,s_y,angle),('skew',s_f,0),('linear',angle,rot))
                model.set_radial_velocity(rvel)
                univ.add_component('observed-'+mol, model)
        """
        Returns a SpectralCube object where all the sources within the FOV and BW are projected.

        This function needs the following parameters:
        - name    : name of the cube
        - alpha   : right-ascension center
        - delta   : declination center
        - freq    : spectral center (frequency)
        - ang_res : angular resolution
        - ang_fov : angular field of view
        - spe_res : spectral resolution
        - spe_bw  : spectral bandwidth
        """
        alpha = 0.0
        delta = 0.0
        #On the Band 9 (602 - 720 Ghz),
        # a sample of 4 Ghz with resolution 1 Mhz:
        # Sample: [602 - 606]
        freq = 604000 + 4000*freq_iteration
        spe_bw= 4000
        spe_res= 1
        cube = univ.gen_cube('observerd', alpha, delta, freq, 10, 40, spe_res, spe_bw)
        cube_noiseless = univ.gen_cube('noiseless', alpha, delta, freq, 10, 40, spe_res, spe_bw, white_noise=False)

        win_len = int(s_f)
        if (win_len % 2) == 0:
          win_len -= 1


        values_noise = cube.get_spectrum(0.0,0.0)
        values = cube.get_spectrum(2,2)
        values_noiseless = cube_noiseless.get_spectrum(2,2)
        plt.plot(values)
        plt.plot(values_noiseless)
        #std = np.std(savgol_filter(values_noise, win_len, 2))
        std = np.std(values_noise)
        plt.plot(np.zeros(len(values_noise)) + std)

        #values_denoised = savgol_filter(values, win_len, 2)
        values_denoised = values

        plt.plot(values_denoised)

        sensibity = .005
        maxtab, mintab = utils.peakdet.peakdet(values_denoised,sensibity)
        while len(maxtab) > 0:
            max_line_temp = max(maxtab[:,1])

            if (max_line_temp>2*std):
              plt.plot(maxtab[:,0], maxtab[:,1], 'r*');
              max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]
              gaussian_fitted = max_line_temp*utils.utils.gaussian(np.arange(0,4000,1), max_line_freq, s_f)
              values_denoised = values_denoised - gaussian_fitted
              # plt.plot(values_denoised)
              maxtab, mintab = utils.peakdet.peakdet(values_denoised,sensibity)
            else:
              break

        plt.show()
        break

        # s_f : parameter to detect lines

        range = str(int((freq - spe_bw/2.0)/1000.0)) + " - " + str(int((freq + spe_bw/2.0)/1000.0))
        univ.save_cube(cube,'observed/' + str(Comb) + ' ' + range + '.fits')

