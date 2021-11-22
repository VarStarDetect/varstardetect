from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import pandas as pd

import lightkurve as lk  # for downloading TESS data products


class StarProcessor:

    OUTLIERS = True
    SIGMA = 3
    EPOCH_0 = 0

    def __init__(self) -> None:
        pass

    def data_download(self, target_number, outliers, sigma, dir):

        TICID, Camera, CCD, Tmag, RA, Dec = np.loadtxt(
            f"{dir}", delimiter=',', unpack=True)

        target = int(TICID[target_number])
        # serching target with the lightkurve
        search_result = lk.search_lightcurve(f'TIC {target}')
        lcf = search_result[1].download()  # choosing datafile

        flux = np.array(lcf.pdcsap_flux)
        flux_err = np.array(lcf.pdcsap_flux_err)
        time = np.array(lcf.time.jd)

        # ---------    NaN VALUES    ---------------------------------------------------------------------

        lol = np.argwhere(np.isnan(flux))
        flux = np.delete(flux, lol)
        time = np.delete(time, lol)
        flux_err = np.delete(flux_err, lol)

        # ---------    OUTLIERS    -----------------------------------------------------------------------

        if outliers == True:

            data = {'time': time, 'flux': flux, 'flux_err': flux_err}
            df = pd.DataFrame(data, columns=['time', 'flux', 'flux_err'])

            highest_allowed = df['flux'].mean() + sigma * df['flux'].std()
            lowest_allowed = df['flux'].mean() - sigma * df['flux'].std()

            max_outliers_index = np.array(np.argwhere(
                flux > highest_allowed).flatten(), dtype='int')
            min_outliers_index = np.array(np.argwhere(
                flux < lowest_allowed).flatten(), dtype='int')

            outliers = np.concatenate((max_outliers_index, min_outliers_index))

            time = np.delete(time, outliers)
            flux = np.delete(flux, outliers)
            flux_err = np.delete(flux_err, outliers)

        return time, flux, flux_err

    def period_search(self, time, flux, flux_err, min_period, max_period):
        ls = LombScargle(time, flux, flux_err)

        frequency, power = ls.autopower(
            minimum_frequency=1./max_period, maximum_frequency=1./min_period)

        period = 1./frequency  # period is the inverse of frequency
        # period with most power / strongest signal
        best_period = 1./frequency[np.argmax(power)]

        uncertainty = best_period ** 2 / (time[len(time) - 1] - time[0])

        return best_period, uncertainty

    def angular_frequency(self, best_period, best_period_err):

        omega = 2 * np.pi / best_period

        uncertainty = 2 * np.pi * np.log(best_period) * best_period_err

        return omega, uncertainty

    def harm_fourier_series(self, degree, omega, T_0):

        def func(t):
            y = [1]

            for i in range(1, degree+1):
                n = np.cos(i * omega * (t - T_0))
                y.append(n)

                n = np.sin(i * omega * (t - T_0))
                y.append(n)
            return y

        return func

    def min_squares(self, time, flux, func):

        lf = len(func(1))
        l = len(time)

        M = np.zeros((l, lf))

        for i in range(0, l):
            for j in range(0, lf):
                M[i, j] = func(time[i])[j]

        A = np.matmul(M.transpose(), M)
        B = np.matmul(M.transpose(), flux.transpose())

        param = np.linalg.solve(A, B)

        return param

    def min_squares_weighted(self, time, flux, flux_err, func):

        flux = np.asarray(flux)  # The flux is transformed into a matrix

        lf = len(func(1))
        l = len(time)

        M = np.zeros((l, lf))

        for i in range(0, l):
            for j in range(0, lf):
                M[i, j] = func(time[i])[j]

        W = np.zeros((l, l))

        for i in range(0, l):
            W[i, i] = (1 / flux_err[i]) ** 2

        A = np.matmul(np.matmul(M.transpose(), W), M)
        B = np.matmul(np.matmul(M.transpose(), W), flux.transpose())

        param = np.linalg.solve(A, B)

        cvm = np.linalg.inv(A)
        err = np.sqrt(abs(np.diag(cvm)))

        return param, err  # ¿Meter excepcion si det(A) == 0?

    def evaluation(self, func, time, flux, flux_err):

        param, param_err = self.min_squares_weighted(
            time, flux, flux_err, func)

        l = len(time)
        lf = len(func(1))

        y = np.zeros((1, l))

        for i in range(0, l):
            y[0, i] = np.dot(param, func(time[i]))

        y = y.flatten()

        y_err = np.zeros((1, l))

        for j in range(0, l):
            for i in range(0, lf):
                y_err[0, j] = y_err[0, j] + \
                    (func(time[j])[i] * param_err[i]) ** 2

        y_err = np.sqrt(y_err)
        y_err = y_err.flatten()

        return y, y_err

    def chi2_reduced(self, func, time, flux, flux_err):

        deviations = flux - self.evaluation(func, time, flux, flux_err)[0]
        chi2 = np.sum((deviations / flux_err) ** 2)
        return chi2/(len(flux) - len(func(1)))

    def amplitude_func(self, fitted, fitted_err):

        amplitude = max(fitted) - min(fitted)

        amplitude_err = fitted_err[np.argwhere(fitted == max(
            fitted))]**2 + fitted_err[np.argwhere(fitted == min(fitted))]**2
        amplitude_err = np.sqrt(amplitude_err)

        return amplitude, float(amplitude_err[0])

    def best_fit(self, time, flux, flux_err, epoch, omega):
        for i in range(1, 20):

            func = self.harm_fourier_series(i, omega, epoch)

            X2 = self.chi2_reduced(func, time, flux, flux_err)
            delta = abs(1 - X2)

            if i == 1:
                minChi = X2
                mindelta = delta
                mins = i

            elif delta <= mindelta:
                minChi = X2
                mindelta = delta
                mins = i

        return minChi, mins
