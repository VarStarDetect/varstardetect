from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

from astropy.timeseries import LombScargle
from scipy.signal import chirp, find_peaks, peak_widths
from IPython.display import display
import pandas as pd

import lightkurve as lk  # for downloading TESS data product
import csv


def data_download(target_number, outliers, sigma, directory):
    """                
                            data_download DOCUMENTATION
    ---------------------------------------------------------------------
    Downloads time, flux and flux uncertainty value for a given star in 
    the previously downloaded sector list using the lighkurve package.
    ---------------------------------------------------------------------
    INPUTS:     - target_number: position of the target in csv file.
                - outliers: boolean, trim outliers.
                - sigma: number of sigma (outliers).
                - directory: directory with TESS sector 1 files csv. You can
                  download from
                  https://tess.mit.edu/observations/target-lists/
            -------------------------------------------------------------
    OUTPUTS:    - time: julian date.
                - flux: flux of star in e-/s.
                - flux_err: flux uncertainty of star in e-/s.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Imports targets csv files.
                2. Searches for target (search_result).
                3. Downloads target data from the the datapackage selected
                   [no. in square brackets]
                4. Deletes NaN values.
                5. Deletes sigma outliers
                6. Creates arrays with the values given.
    ----------------------------------------------------------------------
    """

    TICID, Camera, CCD, Tmag, RA, Dec = np.loadtxt(
        directory, delimiter=',', unpack=True)

    target = int(TICID[target_number])
    # serching target with the lightkurve
    search_result = lk.search_lightcurve(f'TIC {target}')
    lcf = search_result[1].download()  # choosing datafile

    flux = np.array(lcf.flux)
    flux_err = np.array(lcf.flux_err)
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


def period_search(time, flux, flux_err, min_period, max_period):
    """                
                            period_search DOCUMENTATION
    ---------------------------------------------------------------------
    Obtains the period and its uncertinty of a periodic signal of two 
    dimensional data set using the lomb Scargle method. based on 
    astropy.timeseries. For more information visit:
    https://docs.astropy.org/en/stable/api/astropy.timeseries.LombScargle.html
    ---------------------------------------------------------------------
    INPUTS:     - func: defined harmonic series.
                - time: julian date.
                - flux: flux of star in e-/s.
                - flux_err: flux uncertainty of star in e-/s.
                - min_period: period value to start search.
                - max_period: period value to stop search.
            -------------------------------------------------------------
    OUTPUTS:    - period: period.
                - uncertainty: period uncertainty.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Calcualtes the period and its uncertainty using the 
                   astropy.timeseries Lomb Scargle method.
    ----------------------------------------------------------------------
    """
    ls = LombScargle(time, flux, flux_err)

    frequency, power = ls.autopower(
        minimum_frequency=1./max_period, maximum_frequency=1./min_period)

    # period with most power / strongest signal
    best_period = 1./frequency[np.argmax(power)]

    uncertainty = best_period ** 2 / (time[len(time) - 1] - time[0])

    return best_period, uncertainty


def angular_frequency(time, flux, flux_err, min_period, max_period):
    """                
                            angular_frequency DOCUMENTATION
    ---------------------------------------------------------------------
    Obtains the angular frequency of a periodic sinal in a two 
    dimensional data set.
    ---------------------------------------------------------------------
    INPUTS:     - func: defined harmonic series.
                - time: julian date.
                - flux: flux of star in e-/s.
                - flux_err: flux uncertainty of star in e-/s.
                - min_period: period value to start search.
                - max_period: period value to stop search.
            -------------------------------------------------------------
    OUTPUTS:    - omega: angular frequency.
                - uncertainty: angular frequency uncertainty.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Calcualtes the period and its uncertainty.
                2. Calcuales the angular frequency.
                3. Calculates the uncertainty in the angular frequency.
    ----------------------------------------------------------------------
    """

    best_period, best_period_err = period_search(
        time, flux, flux_err, min_period, max_period)

    omega = 2 * np.pi / best_period

    uncertainty = 2 * np.pi * np.log(best_period) * best_period_err

    return omega, uncertainty


def harm_fourier_series(s, omega, T_0):
    """                
                        harm_fourier_series DOCUMENTATION
    ---------------------------------------------------------------------
    Obtains a function which outputs the terms in an array of the 
    s-degree polynomial.
    ---------------------------------------------------------------------
    INPUTS:     - s: degree of the polynomial.
                - omega: angual frequency of the star.
                - T_0: initial epoch of the star.
            -------------------------------------------------------------
    OUTPUTS:    - func: the array-function of the terms of the polynomial.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Defines the function func.
                2. Calculates the first s terms of the polynomial.
                3. Appends them to the output of func.
    ----------------------------------------------------------------------
    """

    def func(t):
        y = [1]

        for i in range(1, s+1):
            n = np.cos(i * omega * (t - T_0))
            y.append(n)

            n = np.sin(i * omega * (t - T_0))
            y.append(n)
        return y

    return func


def min_squares(time, flux, func):
    """                
                            min_squares DOCUMENTATION
    ---------------------------------------------------------------------
    Obtains the fitting parameters for a given function using the unweighted 
    least squares method. 
    ---------------------------------------------------------------------
    INPUTS:     - func: defined harmonic series
                - time: julian date
                - flux: flux of star in e-/s
            -------------------------------------------------------------
    OUTPUTS:    - param: fitting parameters
    ----------------------------------------------------------------------
    PROCEDURE:
                1. 
    ----------------------------------------------------------------------
    """

    lf = len(func(1))
    l = len(time)

    M = np.zeros((l, lf))

    for i in range(0, l):
        for j in range(0, lf):
            M[i, j] = func(time[i])[j]

    A = np.matmul(M.transpose(), M)
    B = np.matmul(M.transpose(), flux.transpose())

    param = np.linalg.solve(A, B)
    #param = np.linalg.lstsq(A, B, rcond=None)[0]

    # IMPORTANTE se hace esto para minimizar A*param = B si el sistema dado es singular(tiene det. 0, con lo que infty sol.)
    # https://stackoverflow.com/questions/64527098/numpy-linalg-linalgerror-singular-matrix-error-when-trying-to-solve

    return param


def min_squares_weighted(time, flux, flux_err, func):
    """                
                        min_squares_weighted DOCUMENTATION
    ---------------------------------------------------------------------
    Obtains the fitting parameters for a given function using the weighted
    least squares method. 
    ---------------------------------------------------------------------
    INPUTS:     - func: defined harmonic series
                - time: julian date
                - flux: flux of star in e-/s
                - flux_err : uncertainty of the flux measurement
            -------------------------------------------------------------
    OUTPUTS:    - param: fitting parameters
                - err: uncertainty of each parameter
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Defines the matrices of the function and weights.
                2. Proposes the normal equations for the weighted
                   least squares fitting method.
                3. Calculates the covariance matrix and the errors of
                   each parameter of the fit. 
    ----------------------------------------------------------------------
    """

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
    #param = np.linalg.lstsq(A, B, rcond=None)[0]

    # IMPORTANTE se hace esto para minimizar A*param = B si el sistema dado es singular(tiene det. 0, con lo que infty sol.)
    # https://stackoverflow.com/questions/64527098/numpy-linalg-linalgerror-singular-matrix-error-when-trying-to-solve

    return param


def evaluation(func, time, flux, flux_err):

    param = min_squares_weighted(time, flux, flux_err, func)

    l = len(time)
    lf = len(func(1))

    y = np.zeros((1, l))  # Se calcula la evaluacion de la funcion

    for i in range(0, l):
        y[0, i] = np.dot(param, func(time[i]))

    y = y.flatten()

    y_err = np.zeros((1, l))  # Se calculan los errores de cada evaluacion

    # Calculo A
    A = np.zeros((lf, lf))

    for a in range(0, lf):
        for b in range(0, lf):

            e = 0
            for k in range(0, l):

                f = func(time[k])
                e = e + f[a]*f[b]    # Esto es A_{aplha,beta}

            A[a, b] = e

    A_inv = np.linalg.inv(A)

    for i in range(0, l):

        sig2 = 0
        for a in range(0, lf):
            for b in range(0, lf):

                ff = func(time[i])
                # Con esto el sum A_{alpha,beta}*f_alpha(t)*f_beta(t)
                sig2 = sig2 + (ff[a]*ff[b])*A_inv[a, b]

        sigstar = 0
        for d in range(1, l):
            sigstar = sigstar + (flux[d] - y[d])**2

        sigstar = sigstar/(l-lf)  # Con esto sigma*^2

        # Con esto se tiene la raiz cuadrada de la formula 4 con grd 0
        y_err[0, i] = np.sqrt(sig2*sigstar)

    y_err = y_err.flatten()

    return y, y_err


# Pa no calcular el error y optimizar chi^2
def evaluation_ini(func, time, flux, flux_err):

    param = min_squares_weighted(time, flux, flux_err, func)

    l = len(time)
    lf = len(func(1))

    y = np.zeros((1, l))  # Se calcula la evaluacion de la funcion

    for i in range(0, l):
        y[0, i] = np.dot(param, func(time[i]))

    y = y.flatten()

    return y


def chi2_reduced(func, time, flux, flux_err):
    """                
                            chi2_reduced DOCUMENTATION
    ---------------------------------------------------------------------
    Calculates the reduced chi^2 parameter for a given harmonic series
    function. 
    ---------------------------------------------------------------------
    INPUTS:     - func: defined harmonic series.
                - time: julian date.
                - flux: flux of star in e-/s .
                - flux_err: flux uncertainty of star in e-/s.
            -------------------------------------------------------------
    OUTPUTS:    - chi2_reduced: reduced chi^2 parameter.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Calculates deviations from the expiramental date and
                   the evaluated values of the fitted function.
                2. Calculates chi^2 parameter.
                3. Calculates degrees of fredom.
                4. Calculates the reduced chi^2 parameter.
    ----------------------------------------------------------------------
    """
    deviations = flux - evaluation_ini(func, time, flux, flux_err)
    chi2 = np.sum((deviations / flux_err) ** 2)
    return chi2/(len(flux) - len(func(1)))


def best_fit(time, flux, flux_err, omega, T_0):
    """                
                            best_fit DOCUMENTATION
    ---------------------------------------------------------------------
    Calculates degree of fourier polynomial that minimises chi^2. Max
    degree 15.
    ---------------------------------------------------------------------
    INPUTS:     - time: julian date.
                - flux: flux of star in e-/s.
                - flux_err: flux uncertainty of star in e-/s.
                - omega: angular frequency of periodic signal.
                - T_0: epoch.
            -------------------------------------------------------------
    OUTPUTS:    - minChi: reduced chi^2 paraemter of the optimal fit.
                - mins: degree of the optimal harmonic polynomial.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Defines harmonic polynomial or every degree (i) in range.
                2. Calculates chi^2 parameter for every func value.
                3. Obtains degree (i) with best chi^2 =~ 1 approximation.
    ----------------------------------------------------------------------
    """

    for i in range(1, 21):

        func = harm_fourier_series(i, omega, T_0)
        X2 = chi2_reduced(func, time, flux, flux_err)
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


def amplitude_func(flux, flux_err, time, mins, omega, T_0):
    """                
                        amplitude_func DOCUMENTATION
    ---------------------------------------------------------------------
    Calculates amplitude light curve of a given star data and its 
    uncertainty.
    ---------------------------------------------------------------------
    INPUTS:     - flux: flux of star in e-/s.
                - time: julian date.
                - mins: degree of harmonic polynomial.
                - omega: angular frequency of periodic signal.
                - T_0: epoch.
            -------------------------------------------------------------
    OUTPUTS:    - amplitude: amplitude of fitted curve.
                - amplitude_err: amplitude uncertainty of fitted curve.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Defines harmonic polynomial.
                2. Fits function to data using least squares and calculates
                   the fitting parameters.
                3. Calculates amplitude
                4. Calculates amplitude uncertainty.
    ----------------------------------------------------------------------
    """

    func = harm_fourier_series(mins, omega, T_0)
    param = min_squares_weighted(time, flux, flux_err, func)
    fitted, fitted_err = evaluation(func, time, flux, flux_err)
    amplitude = max(fitted) - min(fitted)

    amplitude_err = fitted_err[np.argwhere(fitted == max(
        fitted))]**2 + fitted_err[np.argwhere(fitted == min(fitted))]**2
    amplitude_err = np.sqrt(amplitude_err)

    return amplitude, amplitude_err


def light_curve_fit(target_number, outliers, sigma, min_period, max_period, directory):
    """                
                        light_curve_fit DOCUMENTATION
    ---------------------------------------------------------------------
    Function that calculates light curve of a given star data set. Epoch 
    calculated is minimum. 
    ---------------------------------------------------------------------
    INPUTS:     - target_number: position of star on sector csv file.
                - outliers: boolean, trim outliers.
                - sigma: number of sigma (outliers).
                - min_period: min period to start the search.
                - max_period: max period to finish the search.
                - directory: directory with TESS sector 1 files csv. You can
                  download from
                  https://tess.mit.edu/observations/target-lists/
            -------------------------------------------------------------
    OUTPUTS:   Fit values.
                - fitted: 1D numpy array with evaluated flux values.
                - minChi: chi^2 parameter of the least-squared fit.
                - mins: degree of harmonic (fourier) polynomial.

                Periodicity values.
                - period: period of variable star.
                - period_err: period uncertainty.

                Amplitude values.
                - amplitude: value of the amplitude of the fitted function.
                - amplitude_err: amplitude uncertainty.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Downloads data from TESS database with lightkurve.
                2. Calculates star's epoch from minimum.
                3. Applies Lomb-Scargle algorithm to calculate period.
                4. Calculates angluar frequency (2*pi/T).
                5. Applies least squared method on fourier series to find
                   the degree of polynomial that minimises chi^2.
                6. Defines harmonic polynomial of degree mins.
                7. Evaluates time values with the harmonic pol. defined.
                8. Calculates amplitude and its uncertainty.
    ----------------------------------------------------------------------
    """

    time, flux, flux_err = data_download(
        target_number, outliers, sigma, directory)
    T_0 = 0
    period, period_err = period_search(
        time, flux, flux_err, min_period, max_period)
    omega, omega_err = angular_frequency(
        time, flux, flux_err, min_period, max_period)
    minChi, mins = best_fit(time, flux, flux_err, omega, T_0)
    func = harm_fourier_series(mins, omega, T_0)
    fitted, fitted_err = evaluation(func, time, flux, flux_err)
    amplitude, amplitude_err = amplitude_func(
        flux, flux_err, time, mins, omega, T_0)

    return time, flux, flux_err, fitted, minChi, mins, period, period_err, amplitude, amplitude_err


def phase_plot_fit(target_number, outliers, sigma, min_period, max_period, directory):
    """                
                        phase_plot_fit DOCUMENTATION
    ---------------------------------------------------------------------
    Function that calculates phase plot of a given star data set. 
    ---------------------------------------------------------------------
    INPUTS:     - target_number: position of star on sector csv file.
                - outliers: boolean, trim outliers.
                - sigma: number of sigma (outliers).
                - min_period: min period to star the search.
                - max_period: max period to finish the search.
                - directory: directory with TESS sector 1 files csv. You can
                  download from
                  https://tess.mit.edu/observations/target-lists/
            -------------------------------------------------------------
    OUTPUTS:   Data values
                - time
                - phase
                - flux,
                - flux_err

                Fit values.
                - fitted: 1D numpy array with evaluated flux values
                - minChi: chi^2 parameter of the least-squared fit
                - mins: degree of harmonic (fourier) polynomial

                Periodicity values.
                - period: period of variable star
                - period_err: period uncertainty

                Amplitude values.
                - amplitude: value of the amplitude of the fitted function.
                - amplitude_err: uncertainty of the amplitude.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Downloads data from TESS database with lightkurve.
                2. Calculates star's epoch from minimum.
                3. Applies Lomb-Scargle algorithm to calculate period.
                4. Calculates angluar frequency (2*pi/T).
                5. Calculates phase
                6. Applies least squared method on fourier series to find
                   the degree of polynomial that minimises chi^2.
                7. Defines harmonic polynomial of degree mins.
                8. Evaluates time values with the harmonic pol. defined.
                9. Calculates amplitude and its uncertainty.
    ----------------------------------------------------------------------
    """

    time, flux, flux_err = data_download(
        target_number, outliers, sigma, directory)

    T_0 = time[np.argmin(flux)]

    period, period_err = period_search(
        time, flux, flux_err, min_period, max_period)
    omega, omega_err = angular_frequency(
        time, flux, flux_err, min_period, max_period)

    phase = (time - T_0) / period % 1

    minChi, mins = best_fit(phase, flux, flux_err, omega, T_0)

    func = harm_fourier_series(mins, omega, T_0)

    fitted, fitted_err = evaluation(func, phase, flux, flux_err)
    amplitude, amplitude_err = amplitude_func(
        flux, flux_err, phase, mins, omega, T_0)

    return time, phase, flux, flux_err, fitted, minChi, mins, period, period_err, amplitude, amplitude_err


def phase_plot(target_number, outliers, sigma, min_period, max_period, directory):
    """                
                        phase_plot DOCUMENTATION
    ---------------------------------------------------------------------
    Function that calculates light curve of a given star data set. 
    ---------------------------------------------------------------------
    INPUTS:     - target_number: position of star on sector csv file.
                - outliers: boolean, trim outliers.
                - sigma: number of sigma (outliers).
                - min_period: min period to star the search.
                - max_period: max period to finish the search.
                - directory: directory with TESS sector 1 files csv. You can
                  download from
                  https://tess.mit.edu/observations/target-lists/
            -------------------------------------------------------------
    OUTPUTS:   - phase plot 
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Applies phase_plot_fit
                2. Returns plot
    ----------------------------------------------------------------------
    """
    TICID, Camera, CCD, Tmag, RA, Dec = np.loadtxt(
        directory, delimiter=',', unpack=True)
    time, phaseA, flux, flux_err, fitted, minChi, mins, period, period_err, amplitude, amplitude_err = phase_plot_fit(
        target_number, outliers, sigma, min_period, max_period, directory)

    phase = np.sort(np.concatenate((phaseA, phaseA + 1)))
    residuals = flux - fitted
    fig, phase_plot = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={
                                   'height_ratios': [3, 1]})

    fig.suptitle(f"TIC {int(TICID[target_number])}")
    phase_plot[0].plot(phase, np.concatenate((flux, flux)), "+")
    phase_plot[0].plot(phase, np.concatenate((fitted, fitted)))
    phase_plot[1].plot(phase, np.concatenate((residuals, residuals)), "+")
    fig.savefig(f'{int(TICID[target_number])}.png')
    return phase_plot


def amplitude_test(min, max, amp, directory):
    """                
                       amplitude_test DOCUMENTATION
    ---------------------------------------------------------------------
    Detects variable stars with amplitude higher than threshold.
    ---------------------------------------------------------------------
    INPUTS:     - min: lower star search (TESS) range delimiter
                - max: higher star search (TESS) range delimiter
                - amp: amplitude threshold
                - directory: directory with TESS sector 1 files csv. You can
                  download from
                  https://tess.mit.edu/observations/target-lists/
            -------------------------------------------------------------
    OUTPUTS:    - candidates: 1D numpy array with variable candidate 
                  target IDs
                - chis: 1D numpy array with the chi^2 parameter of each
                  approximation.
                - degree: 1D numpy array with the degree of the optimal
                  degree of the approximation.
                - periods: 1D numpy array with the period of each 
                  approimating function.
                - period_errors: 1D numpy array with the period 
                  uncertainty for each candidate.
                - amplitudes: 1D numpy array with the amplitude of each
                  approximation.
                - amplitude_errors: 1D numpy array with the uncertainty 
                  of the amplitude of each candidate.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Calculates amplitude for an observed star.
                2. Calculates if amplitude is bigger than threshold.
                3. Returns candidates and their characteristics.
    ----------------------------------------------------------------------
    """

    TICID, Camera, CCD, Tmag, RA, Dec = np.loadtxt(
        directory, delimiter=',', unpack=True)
    candidates = []
    chis = []
    degrees = []
    periods = []
    period_errors = []
    amplitudes = []
    amplitude_errors = []
    ticid = []

    for i in range(min, max):

        time, flux, flux_err, fitted, minChi, mins, period, period_err, amplitude, amplitude_err = light_curve_fit(
            i, True, 4, 0.1, 27., directory)

        if amplitude >= amp:

            candidates.append(i)
            ticid.append(int(TICID[i]))

            chis.append(minChi)
            degrees.append(mins)
            periods.append(period)
            period_errors.append(period_err)
            amplitudes.append(amplitude)
            amplitude_errors.append(amplitude_err)

    return candidates, ticid, chis, degrees, periods, period_errors, amplitudes, amplitude_errors


def csv_amplitude_test(file, min, max, amp, directory):
    """                
                       amplitude_test DOCUMENTATION
    ---------------------------------------------------------------------
    Detects variable stars with amplitude higher than threshold.
    ---------------------------------------------------------------------
    INPUTS: 
                - file: csv where to store the results.
                - min: lower star search (TESS) range delimiter
                - max: higher star search (TESS) range delimiter
                - amp: amplitude threshold
                - directory: directory with TESS sector 1 files csv. You can
                  download from
                  https://tess.mit.edu/observations/target-lists/
            -------------------------------------------------------------
    OUTPUTS:    - candidates: 1D numpy array with variable candidate 
                  target IDs
                - chis: 1D numpy array with the chi^2 parameter of each
                  approximation.
                - degree: 1D numpy array with the degree of the optimal
                  degree of the approximation.
                - periods: 1D numpy array with the period of each 
                  approimating function.
                - period_errors: 1D numpy array with the period 
                  uncertainty for each candidate.
                - amplitudes: 1D numpy array with the amplitude of each
                  approximation.
                - amplitude_errors: 1D numpy array with the uncertainty 
                  of the amplitude of each candidate.
                - csv file containing the results.
    ----------------------------------------------------------------------
    PROCEDURE:
                1. Calculates amplitude for an observed star.
                2. Calculates if amplitude is bigger than threshold.
                3. Returns candidates and their characteristics.
    ----------------------------------------------------------------------
    """
    candres, ticid, chisres, degres, perres, pererrres, ampres, amperrres = amplitude_test(
        min, max, amp, directory)
    with open(file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',', quotechar='"')
        writer.writerow(['Candidate', 'TICID', 'Chi', 'Degree', 'Period',
                         'Period_error', 'Amplitude', 'Amplitude_error'])
        for i in range(len(candres)):
            writer.writerow([str(candres[i]), str(ticid[i]), str(chisres[i]), str(degres[i]), str(
                perres[i]), str(pererrres[i]), str(ampres[i]), str(amperrres[i])])
    return candres, ticid, chisres, degres, perres, pererrres, ampres, amperrres
