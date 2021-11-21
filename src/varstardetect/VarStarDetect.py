from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import pandas as pd
import csv
import lightkurve as lk  # for downloading TESS data products

from Star import Star
from StarProcessor import StarProcessor
import threading


class VarStarDetect:
    # Constants
    SECTOR_NUMBER = 7
    MIN_PERIOD = 0.1
    MAX_PERIOD = 27
    MAX_THREADS = 5

    # fields

    
    __sector = 1
    
    
    directory = "C:/Users/UO282276/OneDrive - Universidad de Oviedo/Documentos/Programacion/VSC/VarStarDetect/VarStarDetect 1.1.0/src/varstardetect/Targets/sector_1.csv"

    def __init__(self, sector) -> None:
        
        self.__sp = StarProcessor()

        # Selection of sector default 1
        if(sector <= 0 or sector > self.SECTOR_NUMBER):
            self.__sector = 1
        else:
            self.__sector = sector

        self.__processedStars = list()

        pass

    def processStarMT(self, candidate, min, max):
        
        time, flux, flux_err = self.__sp.data_download(
            candidate, StarProcessor.OUTLIERS, StarProcessor.SIGMA, self.directory)

        period, period_err = self.__sp.period_search(
        time, flux, flux_err, min, max)
        
        omega, omega_err = self.__sp.angular_frequency(
            period, period_err)
        minChi, mins = self.__sp.best_fit(time, flux, flux_err, StarProcessor.EPOCH_0, omega)

        func = self.__sp.harm_fourier_series(mins, omega, StarProcessor.EPOCH_0)

        fitted, fitted_err = self.__sp.evaluation(func, time, flux, flux_err)

        amplitude, amplitude_err = self.__sp.amplitude_func(fitted, fitted_err)

        return Star(candidate, minChi, mins, period, period_err, amplitude, amplitude_err)

    def processStar(self, candidate, min, max):
        
        time, flux, flux_err = self.__sp.data_download(
            candidate, StarProcessor.OUTLIERS, StarProcessor.SIGMA, self.directory)

        period, period_err = self.__sp.period_search(
        time, flux, flux_err, min, max)
        
        omega, omega_err = self.__sp.angular_frequency(
            period, period_err)
        minChi, mins = self.__sp.best_fit(time, flux, flux_err, StarProcessor.EPOCH_0, omega)

        func = self.__sp.harm_fourier_series(mins, omega, StarProcessor.EPOCH_0)

        fitted, fitted_err = self.__sp.evaluation(func, time, flux, flux_err)

        amplitude, amplitude_err = self.__sp.amplitude_func(fitted, fitted_err)

        return Star(candidate, minChi, mins, period, period_err, amplitude, amplitude_err)

    def analize(self, min, max, amp):
        for i in range(min, max):
            aux_Star = self.processStar(i, self.MIN_PERIOD, self.MAX_PERIOD)
            if aux_Star.getAmplitude() >= amp:
                self.__processedStars.append(aux_Star)
                print(aux_Star)
            else:
                print("CANDIDATE: " + str(aux_Star.getCandidate()) + " not elegible")

    def analizeMT(self, min, max, amp):
        candidates = (max - min) / self.MAX_THREADS-1
        rest = (max - min) %  self.MAX_THREADS-1

        t1 = threading.Thread(target=self.analize, args=(int(min),int(min+candidates),amp))
        t2 = threading.Thread(target=self.analize, args=(int(min+candidates+1),int(min+(candidates*2)),amp))
        t3 = threading.Thread(target=self.analize, args=(int(min+(candidates*2)+1),int(max-(candidates*3)),amp))
        t4 = threading.Thread(target=self.analize, args=(int(min+(candidates*3)+1),int(min+(candidates*4)),amp))
        t5 = threading.Thread(target=self.analize, args=(int(min+(candidates*4)+1),int(max),amp))

        print(1)
        t1.start()
        print(1)
        t2.start()
        print(1)
        t3.start()
        print(1)
        t4.start()
        print(1)
        t5.start()
        print(1)

        t1.join()
        t2.join()
        t3.join()
        t4.join()
        t5.join()

        self.__processedStars.sort()


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


