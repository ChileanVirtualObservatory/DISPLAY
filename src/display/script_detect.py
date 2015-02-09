import utils.peakdet
import utils.utils
from scipy.signal import savgol_filter


if __name__ != "script_detect":

sigma = utils.fwhm2sigma(freq, s_f)
        win_len = int(sigma)
        if (win_len % 2) == 0:
          win_len -= 1

          data = cube.hdulist[0]
          rows = data.shape[1]
          rows = data.shape[2]

          for xi in xrange(1, rows):
            for yi in xrange(1, cols)

            values_noise = cube.get_spectrum(0.0,0.0)
            values = cube.get_spectrum(xi,yi)
            # teo_freq = teo_cube.get_spectrum()
            # plt.plot(teo_freq*np.max(values), color='b')

            poli_order = 3
            # values_noiseless = cube_noiseless.get_spectrum(xi,yi)
            # plt.plot(values, color='g')
            # plt.plot(values_noiseless, color='r')
            mean = np.mean(savgol_filter(values_noise, win_len, poli_order))
            std = np.std(savgol_filter(values_noise, win_len, poli_order))
            # std = np.std(values_noise)
            # plt.plot(np.zeros(len(values_noise)) + 3*std, color='c')

            values_denoised = savgol_filter(values, win_len, poli_order)
            # values_denoised = values

            # plt.plot(values_denoised, color='m')

            scan_sensibity = 3*std
            maxtab, mintab = utils.peakdet.peakdet(values_denoised,scan_sensibity)

            observed_lines = pd.DataFrame(np.zeros(len(values)))
            while len(maxtab) > 0:
              max_line_temp = max(maxtab[:,1])

              if (max_line_temp>mean + 3*std):

                max_line_freq = maxtab[maxtab[:,1] == max_line_temp][:,0]
                # plt.plot(max_line_freq, max_line_temp, 'r*');
                gaussian_fitted = (max_line_temp)*utils.utils.gaussian(np.arange(0,4000,1), max_line_freq, s_f)
                # values_denoised = values_denoised - gaussian_fitted

                for i in xrange(0, len(values_denoised)):
                  values_denoised[i] = max(values_denoised[i] - gaussian_fitted[i], mean)

                  observed_lines[max_line_freq] = 1

                  # plt.plot(values_denoised, color='y')
                  maxtab, mintab = utils.peakdet.peakdet(values_denoised,scan_sensibity)
                  else:
                    break


                    plt.show()

                    break

                    # s_f : parameter to detect lines
