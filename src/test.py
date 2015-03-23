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
