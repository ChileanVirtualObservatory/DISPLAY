#This file is part of ChiVO, the Chilean Virtual Observatory
#A project sponsored by FONDEF (D11I1060)
#Copyright (C) 2015 Universidad Tecnica Federico Santa Maria Mauricio Solar
#                                                            Marcelo Mendoza
#                   Universidad de Chile                     Diego Mardones
#                   Pontificia Universidad Catolica          Karim Pichara
#                   Universidad de Concepcion                Ricardo Contreras
#                   Universidad de Santiago                  Victor Parada
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

hdu_list = fits.open("observed_cube.fits")
j = 0

data = np.array(hdu_list[0].data[:,1,1])
noise = np.array(hdu_list[0].data[:,0,0])
hdu_list.close()
index = np.arange(cube_params['freq'] - int(cube_params['spe_bw']/2),
                     cube_params['freq'] + int(cube_params['spe_bw']/2),
                     cube_params['spe_res'])
plt.plot(index, data, 'r')
plt.plot(index, np.ones(4000)*(np.mean(noise) + 3*np.std(noise)), 'g')

for f in lines.index:
        for line in lines[f]:
          plt.plot(f, lines[f][line], 'bs')
          j = j + 1
plt.show()
