"""
This file is part of ChiVO
Copyright (C) Andres Riveros

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
from asidopy import *
import matplotlib.pyplot as plt
import random
import math
import sys
from itertools import combinations
import numpy as np


# Without redshift (Rvel = 0)
# Temp 300 Kelvin
rvel = 0.0
temp = 300.0


# Molecules and its respective isotopes.
""
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
""

isolist = ('COv=0','COv=1','13COv=0','C18O','C17O','13C17O','13C18O',
           'N2H+v=0', 'N2D+', '15NNH+', 'N15NH+',
          'CNv=0', '13CN', 'C15N',
          'HCNv=0', 'HCNv2=1', 'HCNv2=2','HCN v3=1', 'HC15Nv=0', 'H13CNv2=1', 'H13CNv=0', 'HCNv1=1',
          'HCNv3=1', 'DCNv=0', 'DCNv2=1', 'HCNv2=4', 'HCNv2=1^1-v2=4^0',
          'CSv=0', '13C34Sv=0', 'C36Sv=0', 'C34Sv=0', 'CSv=1-0',
          '13CSv=0', 'C33Sv=0', 'CSv=1', 'C34Sv=1',
          'CCS', 'C13CS', '13CCS', 'CC34S'
          'H2S', 'H2S', 'H234S', 'D2S',
          'H2CS', 'H213CS', 'H2C34S',
          'SO2v=0', '33SO2', '34SO2v=0', 'SO2v2=1',
          'OS18O', 'OS17O',
          'H2CO', 'H2C18O', 'H213CO',
          'HCO+v=0', 'HC18O+', 'HC17O+', 'H13CO+',
          'HC5Nv=0', 'HC5Nv11=1', 'HCC13CCCN', 'HCCCC13CN',
          'HCCC13CCN', 'H13CCCCCN', 'HC13CCCCN'
          'CH3OHvt=0', '13CH3OHvt=0', 'CH318OH', 'CH3OHvt=1',
          '13CH3OHvt=1')



log = sys.stdout
dbpath = 'ASYDO'

# isolist = set(molist)

# AllPossComb : All possible combinations og molecules from isolist
# AllPossComb = sum(map(lambda r: list(combinations(isolist, r)), range(1, len(isolist)+1)), [])

for freq_iteration in range(0, 1):


    # for Comb in AllPossComb:
    for mol in isolist:

        univ=vu.Universe(log)

        #if len(Comb) != 1:
        #    break

        # print Comb

        # for key in Comb:

            # for mol in molist[key]:
      # for mol in Comb:

        univ.create_source('observed-'+mol, 0.0, 0.0)
        s_x=100
        s_y=70
        rot=80
        s_f=85
        angle=math.pi/2.0
        model=vu.IMCM(log,dbpath,mol,temp,('normal',s_x,s_y,angle),('skew',s_f,0),('linear',angle,rot))
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
        freq = 661000 + 4000*freq_iteration
        spe_bw= 4000
        spe_res= 1
        cube = univ.gen_cube('observerd', alpha, delta, freq, 10, 20, spe_res, spe_bw)

        mean = np.mean(cube.get_spectrum())

        if mean != 0:
          print mol
          # plt.plot(cube.get_spectrum(0.0,0.0))
          # plt.show()

        # range = str(int((freq - spe_bw/2.0)/1000.0)) + " - " + str(int((freq + spe_bw/2.0)/1000.0))
        # univ.save_cube(cube,'observed/' + str(Comb) + ' [' + range + '].fits')
