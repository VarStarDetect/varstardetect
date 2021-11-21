from datetime import date



class Star:
    
    
    

    def __init__(self, candidate, minChi, mins, period, period_err, amplitude, amplitude_err) -> None:
        self.__candidate = candidate
        self.__chi = minChi
        self.__degrees = mins
        self.__period = period
        self.__period_err = period_err
        self.__amplitude = amplitude
        self.__amplitude_error = amplitude_err
        self.__name = str(candidate)
        self.__date = date.today()
        self.__code = self.getCode()

        pass

    def __str__(self) -> str:
        aux = "STAR: "
        aux += self.getCode() + "("
        aux += "Candidate nº: " + str(self.getCandidate())
        aux += " - Chi: {:.5f}".format(self.getChi())
        aux += " - Degrees: "+ str(self.getDegrees())
        aux += " - Period: {:.5f}".format(self.getPeriod())
        aux += " - Period Error: {:.5f}".format(self.getPeriodError())
        aux += " - Amplplitude: {:.5f}".format(self.getAmplitude())
        aux += " - Amplitude Error: {:.5f}".format(self.getAmplitudeError())
        aux += " - Process Date: " + str(self.getDate()) + ")"
        return aux

    def getCandidate(self):
        return self.__candidate

    def getName(self):
        return self.__name

    def setName(self, name):
        self.__name = name

    def getChi(self):
        return self.__chi

    def getDegrees(self):
        return self.__degrees

    def getPeriod(self):
        return self.__period

    def getPeriodError(self):
        return self.__period_err

    def getAmplitude(self):
        return self.__amplitude

    def getAmplitudeError(self):
        return self.__amplitude_error

    def getDate(self):
        return self.__date.strftime("%b-%d-%Y")

    def getCode(self):
        aux = "S01"
        aux += "-" + str(self.getCandidate())
        aux += "-" + self.getName().upper()
        return aux
