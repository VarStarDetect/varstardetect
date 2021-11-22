from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import pandas as pd
import csv
import lightkurve as lk  # for downloading TESS data products

from Star import Star
from StarProcessor import StarProcessor


class VarStarDetect:
    # Constants
    SECTOR_NUMBER = 7
    MIN_PERIOD = 0.1
    MAX_PERIOD = 27
    MAX_THREADS = 5

    # fields

    __sector = 'Targets/sector_' + 1 + '.csv'

    def __init__(self, sector) -> None:

        self.__sp = StarProcessor()

        # Selection of sector default 1
        if(sector <= 0 or sector > self.SECTOR_NUMBER):
            self.__sector = self.setSector(1)
        else:
            self.__sector = self.setSector(sector)

        self.__processedStars = list()

        pass

    def setSector(self, sector):
        self.__sector = 'Targets/sector_' + sector + '.csv'

    def processStarMT(self, candidate, min, max):

        time, flux, flux_err = self.__sp.data_download(
            candidate, StarProcessor.OUTLIERS, StarProcessor.SIGMA, self.directory)

        period, period_err = self.__sp.period_search(
            time, flux, flux_err, min, max)

        omega, omega_err = self.__sp.angular_frequency(
            period, period_err)
        minChi, mins = self.__sp.best_fit(
            time, flux, flux_err, StarProcessor.EPOCH_0, omega)

        func = self.__sp.harm_fourier_series(
            mins, omega, StarProcessor.EPOCH_0)

        fitted, fitted_err = self.__sp.evaluation(func, time, flux, flux_err)

        amplitude, amplitude_err = self.__sp.amplitude_func(fitted, fitted_err)

        return Star(candidate, minChi, mins, period, period_err, amplitude, amplitude_err)

    def processStar(self, candidate, min, max):

        time, flux, flux_err = self.__sp.data_download(
            candidate, StarProcessor.OUTLIERS, StarProcessor.SIGMA, self.getCurrentSector())

        period, period_err = self.__sp.period_search(
            time, flux, flux_err, min, max)

        omega, omega_err = self.__sp.angular_frequency(
            period, period_err)
        minChi, mins = self.__sp.best_fit(
            time, flux, flux_err, StarProcessor.EPOCH_0, omega)

        func = self.__sp.harm_fourier_series(
            mins, omega, StarProcessor.EPOCH_0)

        fitted, fitted_err = self.__sp.evaluation(func, time, flux, flux_err)

        amplitude, amplitude_err = self.__sp.amplitude_func(fitted, fitted_err)

        return Star(candidate, minChi, mins, period, period_err, amplitude, amplitude_err)

    def analize(self, min, max, amp, show):
        if show == None:
            show = False

        for i in range(min, max):
            aux_Star = self.processStar(i, self.MIN_PERIOD, self.MAX_PERIOD)
            if aux_Star.getAmplitude() >= amp:
                self.__processedStars.append(aux_Star)
                if show:
                    print(aux_Star)
            else:
                if show:
                    print("CANDIDATE: " +
                          str(aux_Star.getCandidate()) + " not elegible")

    def saveToCSV(self, file):
        with open(file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"')
            writer.writerow(['Candidate', 'Chi', 'Degree', 'Period',
                            'Period_error', 'Amplitude', 'Amplitude_error'])
            for star in self.__processedStars:
                writer.writerow([str(star.getCandidate()), str(star.getChi()), str(star.getDegree()), str(
                    star.getPeriod()), str(star.getPeriodError()), str(star.getAmplitude()), str(star.getAmplitudeError())])

    def printStars(self):
        for star in self.__processedStars:
            print(star)

    def printStarsFrom(self, min, max):
        if len(self.__processedStars) > 1 and min < len(self.__processedStars) and max <= len(self.__processedStars) and min <= max:
            for i in range(min, max):
                print(self.__processedStars[i])

    def getCurrentSector(self):
        return self.__sector
