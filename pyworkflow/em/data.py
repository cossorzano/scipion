# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (info@kinestat.com)
# *
# * Kinestat Pharma
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

import json
import math

from pyworkflow.object import *
from pyworkflow.utils.path import writeMD5, verifyMD5

from constants import *
import numpy as np

class EMObject(OrderedObject):
    """Base object for all EM classes"""
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        
    def __str__(self):
        return self.getClassName()
    
    def getFiles(self):
        """ Get all filePaths """
        return None


class Matrix(Scalar):
    def __init__(self, **args):
        Scalar.__init__(self, **args)
        self._matrix = np.eye(4)
        
    def _convertValue(self, value):
        """Value should be a str with comman separated values
        or a list.
        """
        self._matrix = np.array(json.loads(value))
            
    def getObjValue(self):
        self._objValue = json.dumps(self._matrix.tolist())
        return self._objValue
    
    def setValue(self, i, j, value):
        self._matrix[i, j] = value
        
    def getMatrix(self):
        """ Return internal numpy matrix. """
        return self._matrix
    
    def setMatrix(self, matrix):
        """ Override internal numpy matrix. """
        self._matrix = matrix
        
    def __str__(self):
        return np.array_str(self._matrix)
    
    def _copy(self, other, copyDict, copyId, level=1, ignoreAttrs=[]):
        """ Override the default behaviour of copy
        to also copy array data.
        """
        self.setMatrix(np.copy(other.getMatrix()))
        self._objValue = other._objValue

class PKPDUnit:
    UNIT_TIME_H = 1
    UNIT_TIME_MIN = 2
    UNIT_TIME_SEC = 3
    UNIT_CONC_g_L= 10
    UNIT_CONC_mg_L= 11
    UNIT_CONC_ug_L= 12
    UNIT_CONC_ng_L= 13
    UNIT_CONC_g_mL= 14
    UNIT_CONC_mg_mL= 15
    UNIT_CONC_ug_mL= 16
    UNIT_CONC_ng_mL= 17
    UNIT_CONC_g_uL= 18
    UNIT_CONC_mg_uL= 19
    UNIT_CONC_ug_uL= 20
    UNIT_CONC_ng_uL= 21
    UNIT_WEIGHT_g= 100
    UNIT_WEIGHT_mg= 101
    UNIT_WEIGHT_ug= 102
    UNIT_WEIGHT_ng= 103

    unitDictionary = {
        UNIT_TIME_H: "h",
        UNIT_TIME_MIN: "min",
        UNIT_TIME_SEC: "s",
        UNIT_CONC_g_L: "g/L",
        UNIT_CONC_mg_L: "mg/L",
        UNIT_CONC_ug_L: "ug/L",
        UNIT_CONC_ng_L: "ng/L",
        UNIT_CONC_g_mL: "g/mL",
        UNIT_CONC_mg_mL: "mg/mL",
        UNIT_CONC_ug_mL: "ug/mL",
        UNIT_CONC_ng_mL: "ng/mL",
        UNIT_CONC_g_uL: "g/uL",
        UNIT_CONC_mg_uL: "mg/uL",
        UNIT_CONC_ug_uL: "ug/uL",
        UNIT_CONC_ng_uL: "ng/uL",
        UNIT_WEIGHT_g: "g",
        UNIT_WEIGHT_mg: "mg",
        UNIT_WEIGHT_ug: "ug",
        UNIT_WEIGHT_ng: "ng"
    }

    def __init__(self,unitString=""):
        self.unit = self._fromString(unitString)

    def isTime(self):
        return self.unit>=1 and self.unit<=9

    def isConcentration(self):
        return self.unit>=10 and self.unit<=99

    def isWeight(self):
        return self.unit>=100 and self.unit<=109

    def _fromString(self, unitString):
        if unitString =="":
            return None
        for _unitKey, _unitString in self.unitDictionary.items():
            if _unitString == unitString:
                return _unitKey
        return None

    def _toString(self):
        if self.unit:
            return self.unitDictionary[self.unit]
        else:
            return ""

class PKPDVariable:
    TYPE_NUMERIC = 1000
    TYPE_TEXT = 1001
    TYPE_NAMES = {TYPE_NUMERIC: 'numeric',
                  TYPE_TEXT: 'text'}

    ROLE_TIME = 1010
    ROLE_MEASUREMENT = 1011
    ROLE_LABEL = 1012
    ROLE_NAMES = {ROLE_TIME: 'time',
                  ROLE_MEASUREMENT: 'measurement',
                  ROLE_LABEL: 'label'}

    def __init__(self):
        self.varName = None
        self.varType = None
        self.displayString = ""
        self.role = None
        self.comment = ""

    def parseTokens(self,tokens):
        # t ; h ; numeric[%f] ; time ;

        # Get name
        self.varName = tokens[0].strip()

        # Get units
        unitString = tokens[1].strip()
        self.units = PKPDUnit(unitString)
        if not self.units.unit:
            raise Exception("Unrecognized unit: %s"%unitString)

        # Get type and display
        typeString = tokens[2].strip()
        self.displayString = ""
        leftBracket=typeString.find('[')
        rightBracket=typeString.find(']')
        if leftBracket!=-1 and rightBracket!=-1:
            self.displayString = typeString[leftBracket+1:rightBracket]
            typeString = typeString[0:leftBracket]
        typeString = typeString.lower()
        if typeString=="numeric":
            self.varType = PKPDVariable.TYPE_NUMERIC
            if self.displayString=="":
                self.displayString="%f"
        elif typeString=="text":
            self.varType = PKPDVariable.TYPE_TEXT
            if self.displayString=="":
                self.displayString="%s"
        else:
            raise Exception("Unrecognized type: %s"%typeString)

        # Get role
        roleString = tokens[3].strip().lower()
        if roleString=="time":
            self.role = PKPDVariable.ROLE_TIME
        elif roleString=="measurement":
            self.role = PKPDVariable.ROLE_MEASUREMENT
        elif roleString=="label":
            self.role = PKPDVariable.ROLE_LABEL
        else:
            raise Exception("Unrecognized role: %s"%roleString)

        # Get comment
        self.comment = tokens[4].strip()

    def _printToStream(self,fh):
        displayString = self.displayString.replace("%%","%%%%")

        fh.write("%s ; %s ; %s[%s] ; %s ; %s\n" % (self.varName,
                                                   self.getUnitsString(),
                                                   self.getTypeString(),
                                                   displayString,
                                                   self.getRoleString(),
                                                   self.comment))

    def getTypeString(self):
        return self.TYPE_NAMES.get(self.varType, '')

    def getRoleString(self):
        return self.ROLE_NAMES.get(self.role, '')

    def getUnitsString(self):
        return self.units._toString()


class PKPDDose:
    TYPE_BOLUS = 1
    TYPE_INFUSION = 2

    def __init__(self):
        self.varName = None
        self.doseType = None
        self.doseAmount = None
        self.t0 = None
        self.tF = None
        self.units = None
        self.normalization = ""

    def parseTokens(self,tokens):
        # Dose1; bolus t=0 d=60; mg; dose*weight/1000
        # Dose1; infusion t=0...59 d=1; mg; dose*weight/1000

        # Get name
        self.varName = tokens[0].strip()

        # Get type
        doseString = tokens[1].strip()
        doseTokens = doseString.split(' ')
        if len(doseTokens)!=3:
            raise Exception("Unrecognized dose type %s"%doseString)
        doseTypeString = doseTokens[0].strip().lower()
        timeString = doseTokens[1].strip().lower().split("=")[1]
        if doseTypeString=="bolus":
            self.doseType = PKPDDose.TYPE_BOLUS
            self.t0 = float(timeString)
        elif doseTypeString=="infusion":
            self.doseType = PKPDDose.TYPE_INFUSION
            timeTokens = timeString.split("...")
            self.t0 = float(timeTokens[0])
            self.tF = float(timeTokens[1])
        else:
            raise Exception("Unrecognized dose type %s"%doseTypeString)
        self.doseAmount = doseTokens[2].strip().lower().split("=")[1]

        # Get units
        unitString = tokens[2].strip()
        self.units = PKPDUnit(unitString)
        if not self.units.unit:
            raise Exception("Unrecognized unit: %s"%unitString)
        if not self.units.isWeight():
            raise Exception("After normalization, the dose must be a weight")

    def _printToStream(self,fh):
        fh.write("%s ; %s d=%s; %s\n" % (self.varName,
                                              self.getDoseString(),
                                              self.doseAmount,
                                              self.getUnitsString()))

    def getDoseString(self):
        if self.doseType == PKPDDose.TYPE_BOLUS:
            doseString = "bolus t=%f" % self.t0
        elif self.doseType == PKPDDose.TYPE_INFUSION:
            doseString = "infusion t=%f...%f" % (self.t0, self.tF)
        else:
            doseString = ""

        return doseString

    def getUnitsString(self):
        return self.units._toString()


class PKPDSample:
    def __init__(self):
        self.varName = ""
        self.variableDictPtr = None
        self.doseDictPtr = None
        self.doseList = []
        self.descriptors = None
        self.measurementPattern = None

    def parseTokens(self,tokens,variableDict,doseDict):
        # FemaleRat1; dose=Dose1; weight=207

        # Keep a pointer to variableDict and doseDict
        self.variableDictPtr = variableDict
        self.doseDictPtr = doseDict

        # Get name
        self.varName = tokens[0].strip()

        # Get doses
        doseList = tokens[1].split('=')[1].strip()
        for doseName in doseList.split(','):
            doseName=doseName.strip()
            if doseName in doseDict:
                self.doseList.append(doseName)
            else:
                raise Exception("Unrecognized dose %s"%doseName)

        # Get rest of variables
        self.descriptors = {}
        for n in range(2,len(tokens)):
            if '=' in tokens[n]:
                varTokens = tokens[n].split('=')
                varName  = varTokens[0].strip()
                varValue = varTokens[1].strip()
                if varName in variableDict:
                    varPtr = variableDict[varName]
                    if varPtr.role != PKPDVariable.ROLE_LABEL:
                        raise Exception("Samples can only use role variables")
                    self.descriptors[varName] = varValue

        self.measurementPattern = []

    def addMeasurementPattern(self,tokens):
        self.measurementPattern = []
        for n in range(1,len(tokens)):
            varName = tokens[n].strip()
            if varName in self.variableDictPtr:
                self.measurementPattern.append(varName)
                setattr(self,"measurement_%s"%varName,[])
            else:
                raise Exception("Unrecognized variable %s"%varName)

    def addMeasurement(self,line):
        tokens = line.split()
        if len(tokens)!=len(self.measurementPattern):
            raise Exception("Not enough values to fill measurement pattern")
        for n in range(0,len(tokens)):
            ok=True
            varName = self.measurementPattern[n]
            if tokens[n]=="NA":
                ok = (self.variableDictPtr[varName].role != PKPDVariable.ROLE_TIME)
            if ok:
                exec("self.measurement_%s.append('%s')"%(varName,tokens[n]))
            else:
                raise Exception("Time measurements cannot be NA")

    def getNumberOfVariables(self):
        return len(self.measurementPattern)

    def getNumberOfMeasurements(self):
        return len(getattr(self,"measurement_%s"%self.measurementPattern[0]))

    def _printToStream(self,fh):
        descriptorString = ""
        for key, value in self.descriptors.iteritems():
            descriptorString +="; %s=%s"%(key,value)
        fh.write("%s; dose=%s %s\n"%(self.varName,",".join(self.doseList),descriptorString))

    def _printMeasurements(self,fh):
        patternString = ""
        for n in range(0,len(self.measurementPattern)):
            patternString += "; %s"%self.measurementPattern[n]
        fh.write("%s %s\n"%(self.varName,patternString))
        if len(self.measurementPattern)>0:
            aux=getattr(self,"measurement_%s"%self.measurementPattern[0])
            N = len(aux)
            for i in range(0,self.getNumberOfMeasurements()):
                lineString = ""
                for n in range(0,len(self.measurementPattern)):
                    aux=getattr(self,"measurement_%s"%self.measurementPattern[n])
                    lineString += aux[i]+" "
                fh.write("%s\n"%lineString)
        fh.write("\n")

    def getRange(self, varName):
        if varName not in self.measurementPattern:
            return [None, None]
        else:
            aux = getattr(self,"measurement_%s"%varName)
            aux = [x for x in aux if x != "NA"]
            x = np.array(aux, dtype='|S4')
            y = x.astype(np.float)
            return [y.min(),y.max()]

    def getValues(self, varName):
        if varName not in self.measurementPattern:
            return None
        else:
            return getattr(self,"measurement_%s"%varName)

    def getXYValues(self,varNameX,varNameY):
        xl = []
        yl = []
        xs = self.getValues(varNameX)
        ys = self.getValues(varNameY)
        for x, y in izip(xs, ys):
            if x != "NA" and y != "NA":
                xl.append(float(x))
                yl.append(float(y))
        return xl, yl

    def getSampleMeasurements(self):
        return [PKPDSampleMeasurement(self,n) for n in range(0,self.getNumberOfMeasurements())]


class PKPDSampleMeasurement():
    def __init__(self, sample, n):
        self.sample = sample
        self.n = n

    def getValues(self):
        values = []
        for i in range(0,self.sample.getNumberOfVariables()):
            aux=getattr(self.sample,"measurement_%s"%self.sample.measurementPattern[i])
            values.append(aux[self.n])
        return values


class PKPDExperiment(EMObject):
    READING_GENERAL = 1
    READING_VARIABLES = 2
    READING_DOSES = 3
    READING_SAMPLES = 4
    READING_MEASUREMENTS = 5
    READING_A_MEASUREMENT = 6

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnPKPD = String()
        self.general = {}
        self.variables = {}
        self.samples = {}
        self.doses = {}

    def load(self, fnExperiment="", verifyIntegrity=True):
        if fnExperiment!="":
            self.fnPKPD.set(fnExperiment)
        if verifyIntegrity and not verifyMD5(self.fnPKPD.get()):
            raise Exception("The file %s has been modified since its creation"%self.fnPKPD.get())
        fh=open(self.fnPKPD.get(),'r')
        if not fh:
            raise Exception("Cannot open the file "+self.fnPKPD)

        for line in fh.readlines():
            line=line.strip()
            if line=="":
                if state==PKPDExperiment.READING_A_MEASUREMENT:
                    state=PKPDExperiment.READING_MEASUREMENTS
                continue
            if line[0]=='[':
                section = line.split('=')[0].strip().lower()
                if section=="[experiment]":
                    state=PKPDExperiment.READING_GENERAL
                elif section=="[variables]":
                    state=PKPDExperiment.READING_VARIABLES
                elif section=="[doses]":
                    state=PKPDExperiment.READING_DOSES
                elif section=="[samples]":
                    state=PKPDExperiment.READING_SAMPLES
                elif section=="[measurements]":
                    state=PKPDExperiment.READING_MEASUREMENTS
                else:
                    print("Skipping: ",line)

            elif state==PKPDExperiment.READING_GENERAL:
                tokens = line.split('=')
                lhs = tokens[0].strip().lower()
                rhs = tokens[1].strip()
                self.general[lhs]=rhs

            elif state==PKPDExperiment.READING_VARIABLES:
                tokens = line.split(';')
                if len(tokens)!=5:
                    print("Skipping variable: ",line)
                    continue
                varname = tokens[0].strip()
                self.variables[varname] = PKPDVariable()
                self.variables[varname].parseTokens(tokens)

            elif state==PKPDExperiment.READING_DOSES:
                tokens = line.split(';')
                if len(tokens)!=3:
                    print("Skipping dose: ",line)
                    continue
                dosename = tokens[0].strip()
                self.doses[dosename] = PKPDDose()
                self.doses[dosename].parseTokens(tokens)

            elif state==PKPDExperiment.READING_SAMPLES:
                tokens = line.split(';')
                if len(tokens)<2:
                    print("Skipping sample: ",line)
                    continue
                samplename = tokens[0].strip()
                self.samples[samplename] = PKPDSample()
                self.samples[samplename].parseTokens(tokens,self.variables, self.doses)

            elif state==PKPDExperiment.READING_MEASUREMENTS:
                tokens = line.split(';')
                if len(tokens)<3:
                    print("Skipping measurement: ",line)
                    continue
                samplename = tokens[0].strip()
                if samplename in self.samples:
                    self.samples[samplename].addMeasurementPattern(tokens)
                    state=PKPDExperiment.READING_A_MEASUREMENT
                else:
                    print("Skipping measurement: %s"%line)
            elif state==PKPDExperiment.READING_A_MEASUREMENT:
                self.samples[samplename].addMeasurement(line)

        fh.close()

    def write(self, fnExperiment):
        fh=open(fnExperiment,'w')
        self._printToStream(fh)
        fh.close()
        self.fnPKPD.set(fnExperiment)
        writeMD5(fnExperiment)

    def _printToStream(self,fh):
        fh.write("[EXPERIMENT] ===========================\n")
        for key, value in self.general.iteritems():
            fh.write("%s = %s\n"%(key,value))
        fh.write("\n")

        fh.write("[VARIABLES] ============================\n")
        for key, value in self.variables.iteritems():
            value._printToStream(fh)
        fh.write("\n")

        fh.write("[DOSES] ================================\n")
        for key, value in self.doses.iteritems():
            value._printToStream(fh)
        fh.write("\n")

        fh.write("[SAMPLES] ================================\n")
        for key, value in self.samples.iteritems():
            value._printToStream(fh)
        fh.write("\n")

        fh.write("[MEASUREMENTS] ===========================\n")
        for key, value in self.samples.iteritems():
            value._printMeasurements(fh)
        fh.write("\n")

    def getRange(self,varName):
        vmin = None
        vmax = None
        for key, value in self.samples.iteritems():
            vmini, vmaxi = value.getRange(varName)
            if vmin==None or vmini<vmin:
                vmin = vmini
            if vmax==None or vmaxi>vmax:
                vmax = vmaxi
        return [vmin,vmax]

class PKPDModel:
    def __init__(self):
        self.fnExperiment = None
        self.parameters = None

    def setExperiment(self, experiment):
        self.experiment = experiment
        self.fnExperiment = experiment.fnPKPD

    def setXVar(self, x):
        if not x in self.experiment.variables:
            raise Exception("Cannot find %s as a variable in the experiment"%x)
        self.xName = x
        self.xRange = self.experiment.getRange(x)

    def setYVar(self, y):
        if not y in self.experiment.variables:
            raise Exception("Cannot find %s as a variable in the experiment"%x)
        self.yName = y
        self.yRange = self.experiment.getRange(y)

    def setXYValues(self, x, y):
        self.x = np.array(x)
        self.y = np.array(y)
        self.ylog = np.log(self.y)

    def forwardModel(self, parameters, x=None):
        pass

    def getNumberOfParameters(self):
        pass

    def prepare(self):
        pass

    def printSetup(self):
        pass

    def getEquation(self):
        pass

    def getDescription(self):
        pass

    def setParameters(self, parameters):
        self.parameters = parameters

    def areParametersSignificant(self, lowerBound, upperBound):
        """
        :param lowerBound and upperBound: a numpy array of parameters
        :return: a list of string with "True", "False", "NA", "Suspicious"
        """
        pass

    def setConfidenceInterval(self,lowerBound,upperBound):
        pass


class PKPDExponentialModel(PKPDModel):
    def forwardModel(self, parameters, x=None):
        if x==None:
            x=self.x
        self.yPredicted = np.zeros(x.shape[0])
        proceed=True
        for k in range(1,self.Nexp):
            ck = parameters[2*k]
            ck_1 = parameters[2*(k-1)]
            if ck_1<ck:
                proceed=False
                self.yPredicted = -1000*np.ones(x.shape[0])
                break
        if proceed:
            for k in range(0,self.Nexp):
                ck = parameters[2*k]
                lk = parameters[2*k+1]
                self.yPredicted += ck*np.exp(-lk*x)
        return self.yPredicted

    def setNexp(self, Nexp):
        self.Nexp = Nexp

    def getNumberOfParameters(self):
        return 2*self.Nexp

    def getDescription(self):
        return "Sum of exponentials (%s)"%self.__class__.__name__

    def parseBounds(self,inputString,msg):
        if ';' in inputString:
            tokens=inputString.split(';')
            if len(tokens)!=self.Nexp:
                raise Exception("The number of intervals for %s does not match the number of exponential terms"%msg)
        else:
            tokens = [inputString]*self.Nexp
        return eval("["+",".join(tokens)+"]")

    def setBounds(self, cBounds, lambdaBounds):
        if lambdaBounds=="" or lambdaBounds is None:
            self.lambdaBounds = None
        else:
            self.lambdaBounds = self.parseBounds(lambdaBounds,"lambda")

        if cBounds=="" or cBounds is None:
            self.cBounds = None
        else:
            self.cBounds = self.parseBounds(cBounds,"amplitude")

    def prepare(self):
        if self.cBounds == None and self.lambdaBounds == None:
            print("First estimate of 1 exponential term: ")
            p = np.polyfit(self.x,self.ylog,1)
            self.cBounds=[(math.exp(p[1])*0.01,math.exp(p[1])*100.0)]*self.Nexp
            self.lambdaBounds=[(-p[0]*0.01,-p[0]*100.0)]*self.Nexp
            print("Y=%f*exp(-%f*X)"%(math.exp(p[1]),-p[0]))
        elif self.cBounds == None or self.lambdaBounds == None:
            raise Exception("Either give c and lambda bounds or none of them")
        print("Amplitude (c) bounds: "+str(self.cBounds))
        print("Lambda bounds: "+str(self.lambdaBounds))

    def getBounds(self):
        bounds = []
        for i in range(0,self.Nexp):
            bounds.append(self.cBounds[i])
            bounds.append(self.lambdaBounds[i])
        return bounds

    def printSetup(self):
        print("Model: y=sum_i c_i*exp(-lambda_i * x)")
        print("Number of exponentials: "+str(self.Nexp))

    def getEquation(self):
        toPrint="Y="
        for i in range(self.Nexp):
            toPrint+= "+[%f*exp(-%f*X)]"%(self.parameters[2*i],self.parameters[2*i+1])
        return toPrint

    def areParametersSignificant(self, lowerBound, upperBound):
        retval=[]
        for i in range(self.Nexp):
            cLower = lowerBound[2*i]
            cUpper = lowerBound[2*i+1]
            if cLower<0 and cUpper>0:
                retval.append("False")
            elif cLower>0:
                retval.append("True")
            elif cUpper<0:
                retval.append("Suspicious, this term may be negative")
            else:
                retval.append("NA")

            decayLower = lowerBound[2*i+1]
            decayUpper = lowerBound[2*i+1]
            if decayLower<0 and decayUpper>0:
                retval.append("Suspicious, looks like a constant")
            elif decayLower<0:
                retval.append("Suspicious, this term may be unstable")
            elif decayLower>0:
                retval.append("True")
            else:
                retval.append("NA")
        return retval

    def setConfidenceInterval(self,lowerBound,upperBound):
        yPredictedBackup = self.yPredicted.copy()
        self.yPredictedUpper = self.yPredicted.copy()
        self.yPredictedLower = self.yPredicted.copy()
        for i in range(0,int(math.pow(2,2*self.Nexp))):
            pattern = ("{0:0%db}"%(2*self.Nexp)).format(i)
            p = np.where(np.array(list(pattern))=="1",upperBound,lowerBound)
            p = p.astype(np.float)
            if np.sum(p<0)!=0:
                continue
            y = self.forwardModel(p)
            for n in range(len(y)):
                if y[n]<self.yPredictedLower[n]:
                    if y[n]<0:
                        self.yPredictedLower[n]=0
                    else:
                        self.yPredictedLower[n]=y[n]
                if y[n]>self.yPredictedUpper[n]:
                    self.yPredictedUpper[n]=y[n]
        self.yPredicted = yPredictedBackup.copy()

class PKPDOptimizer:
    def __init__(self,model,fitType,goalFunction="RMSE"):
        self.model = model

        self.yTarget = np.array(model.y, dtype=np.float32)
        self.yTargetLogs = np.log10(self.yTarget)
        if fitType=="linear":
            self.takeYLogs = False
            self.takeRelative = False
        elif fitType=="log":
            self.yTarget = self.yTargetLogs
            self.takeYLogs = True
            self.takeRelative = False
        elif fitType=="relative":
            self.takeYLogs = False
            self.takeRelative = True

        if goalFunction=="RMSE":
            self.goalFunction = self.goalRMSE
        else:
            raise Exception("Unknown goal function")

    def getResiduals(self,parameters):
        yPredicted = np.array(self.model.forwardModel(parameters),dtype=np.float32)
        if self.takeYLogs:
            idx = yPredicted>=1e-20
            nonIdx = yPredicted<1e-20
            yPredicted[idx] = np.log10(yPredicted[idx])
            yPredicted[nonIdx] = -100
        diff = self.yTarget - yPredicted
        if self.takeRelative:
            diff = diff/self.yTarget
        return diff

    def goalRMSE(self,parameters):
        rmse = math.sqrt((self.getResiduals(parameters) ** 2).mean())
        return rmse

    def _printFitting(self, x, y, yp):
        for n in range(0,x.shape[0]):
            print("%f %f %f %f"%(x[n],y[n],yp[n],y[n]-yp[n]))
        e = y-yp
        print("------------------------")
        print("Mean error = %f"%np.mean(e))
        print("Std error = %f"%np.std(e))
        self.R2 = (1-np.var(e)/np.var(y))
        print("R2 = %f"%self.R2)
        print("------------------------")

    def printFitting(self):
        yPredicted = np.array(self.model.forwardModel(self.model.parameters),dtype=np.float32)
        print("==========================================")
        print("X     Y    Ypredicted  Error=Y-Ypredicted ")
        print("==========================================")
        self._printFitting(self.model.x, self.model.y, yPredicted)
        print("==================================================================")
        print("X    log10(Y)  log10(Ypredicted)  Error=log10(Y)-log10(Ypredicted)")
        print("==================================================================")
        self._printFitting(self.model.x, np.log10(self.model.y), np.log10(yPredicted))

class PKPDDEOptimizer(PKPDOptimizer):
    def optimize(self):
        from scipy.optimize import differential_evolution
        print("Optimizing with Differential Evolution (DE), a global optimizer")
        self.optimum = differential_evolution(self.goalFunction, self.model.getBounds())
        print("Best DE function value: "+str(self.optimum.fun))
        print("Best DE parameters: "+str(self.optimum.x))
        self.model.setParameters(self.optimum.x)
        print(self.model.getEquation())
        self.printFitting()
        print(" ")
        return self.optimum

class PKPDLSOptimizer(PKPDOptimizer):
    def optimize(self):
        from scipy.optimize import leastsq
        print("Optimizing with Least Squares (LS), a local optimizer")
        print("Initial parameters: "+str(self.model.parameters))
        self.optimum, self.cov_x, self.info, _, _ = leastsq(self.getResiduals, self.model.parameters, full_output=True)
        print("Best LS function value: "+str(self.goalFunction(self.optimum)))
        print("Best LS parameters: "+str(self.optimum))
        self.cov_x *= np.var(self.info["fvec"])

        print("Covariance matrix:")
        print(np.array_str(self.cov_x,max_line_width=120))
        self.model.setParameters(self.optimum)
        print(self.model.getEquation())
        self.printFitting()
        print(" ")
        return self.optimum

    def setConfidenceInterval(self,confidenceInterval):
        from scipy.stats import norm
        nstd = norm.ppf(1-(1-confidenceInterval/100)/2)
        perr = np.sqrt(np.diag(self.cov_x))
        self.lowerBound = self.optimum-nstd*perr
        self.upperBound = self.optimum+nstd*perr

        self.significance = self.model.areParametersSignificant(self.lowerBound,self.upperBound)
        print("Confidence intervals %f%% --------------------------"%confidenceInterval)
        print("ParameterValue   ParameterConfidenceInterval  IsStatisticallySignificant")
        for n in range(0,len(self.optimum)):
            print("%f [%f,%f] %s"%(self.optimum[n],self.lowerBound[n],self.upperBound[n],self.significance[n]))

        self.model.setConfidenceInterval(self.lowerBound,self.upperBound)

class PKPDSampleFit:
    def __init__(self):
        self.sampleName = ""
        self.x = None
        self.y = None
        self.yp = None
        self.yl = None
        self.yu = None
        self.modelEquation = ""
        self.R2 = 0
        self.parameters = None
        self.lowerBound = None
        self.upperBound = None
        self.significance = None

    def _printToStream(self,fh):
        if self.significance == None:
            self.significance = ["Undetermined"]*len(self.parameters)
        fh.write("Sample name: %s\n"%self.sampleName)
        fh.write("Model: %s\n"%self.modelEquation)
        fh.write("R2: %f\n"%self.R2)
        fh.write("Parameter lowerBound upperBound IsStatisticallySignificant -------\n")
        for parameter, lower, upper, significance in izip(self.parameters,self.lowerBound,self.upperBound,\
                                                          self.significance):
            fh.write("%f [%f,%f] %s\n"%(parameter,lower,upper,significance))
        fh.write("X   Y   Ypredicted [Ylower,Yupper] -------\n")
        for x,y,yp,yl,yu in izip(self.x,self.y,self.yp,self.yl,self.yu):
            fh.write("%f %f %f [%f,%f]\n"%(x,y,yp,yl,yu))
        fh.write("\n")

class PKPDFitting(EMObject):
    READING_FITTING_EXPERIMENT = 1
    READING_FITTING_PREDICTOR = 2
    READING_FITTING_PREDICTED = 3
    READING_FITTING_MODEL = 4
    READING_POPULATION = 5
    READING_SAMPLEFITTINGS_NAME = 6
    READING_SAMPLEFITTINGS_MODELEQ = 7
    READING_SAMPLEFITTINGS_R2 = 8
    READING_SAMPLEFITTINGS_PARAMETER_BOUNDS = 9
    READING_SAMPLEFITTINGS_SAMPLE_VALUES = 10

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnFitting = String()
        self.fnExperiment = String()
        self.predictor = None
        self.predicted = None
        self.modelDescription = ""
        self.sampleFits = []

    def write(self, fnFitting):
        fh=open(fnFitting,'w')
        self._printToStream(fh)
        fh.close()
        self.fnFitting.set(fnFitting)
        writeMD5(fnFitting)

    def _printToStream(self,fh):
        fh.write("[FITTING] ===========================\n")
        fh.write("Experiment: %s\n"%self.fnExperiment.get())
        fh.write("Predictor (X): ")
        self.predictor._printToStream(fh)
        fh.write("Predicted (Y): ")
        self.predicted._printToStream(fh)
        fh.write("Model: %s\n"%self.modelDescription)
        fh.write("\n")

        fh.write("[POPULATION PARAMETERS] =============\n")
        i=0
        for sampleFitting in self.sampleFits:
            outputStr = ""
            j=0
            for parameter in sampleFitting.parameters:
                outputStr += "%f "%parameter
                if i==0 and j==0:
                    observations = np.zeros([len(self.sampleFits),len(sampleFitting.parameters)])
                observations[i,j]=parameter
                j+=1
            fh.write(outputStr+"\n")
            i+=1
        fh.write("\n")

        mu=np.mean(observations,axis=0)
        C=np.cov(np.transpose(observations))
        sigma = np.sqrt(np.diag(C))
        fh.write("Mean   parameters  = %s\n"%np.array_str(mu))
        fh.write("Median parameters  = %s\n"%np.array_str(np.median(observations,axis=0)))
        fh.write("Lower bound (95%%, independent Gaussians) = %s\n"%np.array_str(mu-1.96*sigma))
        fh.write("Upper bound (95%%, independent Gaussians) = %s\n"%np.array_str(mu+1.96*sigma))
        fh.write("Covariance matrix  =\n%s\n"%np.array_str(C,max_line_width=120))
        fh.write("\n")

        fh.write("[SAMPLE FITTINGS] ===================\n")
        for sampleFitting in self.sampleFits:
            sampleFitting._printToStream(fh)

    def load(self,fnFitting):
        fh=open(fnFitting)
        if not fh:
            raise Exception("Cannot open %s"%fnFitting)
        if not verifyMD5(fnFitting):
            raise Exception("The file %s has been modified since its creation"%fnFitting)
        self.fnFitting.set(fnFitting)

        for line in fh.readlines():
            line=line.strip()
            if line=="":
                if state==PKPDFitting.READING_SAMPLEFITTINGS_SAMPLE_VALUES:
                     state=PKPDFitting.READING_SAMPLEFITTINGS_NAME
                continue
            if line[0]=='[':
                section = line.split('=')[0].strip().lower()
                if section=="[fitting]":
                    state=PKPDFitting.READING_FITTING_EXPERIMENT
                elif section=="[population parameters]":
                    state=PKPDFitting.READING_POPULATION
                elif section=="[sample fittings]":
                    state=PKPDFitting.READING_SAMPLEFITTINGS_NAME
                else:
                    print("Skipping: ",line)

            elif state==PKPDFitting.READING_FITTING_EXPERIMENT:
                tokens = line.split(':')
                self.fnExperiment.set(tokens[1].strip())
                state = PKPDFitting.READING_FITTING_PREDICTOR

            elif state==PKPDFitting.READING_FITTING_PREDICTOR:
                tokens = line.split(':')
                self.predictor = PKPDVariable()
                self.predictor.parseTokens(tokens[1].split(';'))
                state = PKPDFitting.READING_FITTING_PREDICTED

            elif state==PKPDFitting.READING_FITTING_PREDICTED:
                tokens = line.split(':')
                self.predicted = PKPDVariable()
                self.predicted.parseTokens(tokens[1].split(';'))
                state = PKPDFitting.READING_FITTING_MODEL

            elif state==PKPDFitting.READING_FITTING_MODEL:
                tokens = line.split(':')
                self.modelDescription = tokens[1].strip()
                state = PKPDFitting.READING_POPULATION

            elif state==PKPDFitting.READING_POPULATION:
                continue

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_NAME:
                tokens = line.split(':')
                self.sampleFits.append(PKPDSampleFit())
                self.sampleFits[-1].sampleName = tokens[1].strip()
                state = PKPDFitting.READING_SAMPLEFITTINGS_MODELEQ

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_MODELEQ:
                print("Model line=====",line)
                tokens = line.split(':')
                self.sampleFits[-1].modelEquation = tokens[1].strip()
                state = PKPDFitting.READING_SAMPLEFITTINGS_R2

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_R2:
                tokens = line.split(':')
                self.sampleFits[-1].R2 = float(tokens[1])
                state = PKPDFitting.READING_SAMPLEFITTINGS_PARAMETER_BOUNDS

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_PARAMETER_BOUNDS:
                if line.startswith("Parameter lowerBound upperBound"):
                    self.sampleFits[-1].parameters=[]
                    self.sampleFits[-1].lowerBound=[]
                    self.sampleFits[-1].upperBound=[]
                    self.sampleFits[-1].significance=[]
                elif line.startswith("X   Y   Ypredicted"):
                    state=PKPDFitting.READING_SAMPLEFITTINGS_SAMPLE_VALUES
                    self.sampleFits[-1].x=[]
                    self.sampleFits[-1].y=[]
                    self.sampleFits[-1].yp=[]
                    self.sampleFits[-1].yl=[]
                    self.sampleFits[-1].yu=[]
                else:
                    tokens=line.split()
                    self.sampleFits[-1].parameters.append(float(tokens[0]))
                    self.sampleFits[-1].significance.append(tokens[2])
                    tokens=(tokens[1])[1:-1].split(',')
                    self.sampleFits[-1].lowerBound.append(float(tokens[0]))
                    self.sampleFits[-1].upperBound.append(float(tokens[1]))

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_SAMPLE_VALUES:
                tokens=line.split()
                self.sampleFits[-1].x.append(float(tokens[0]))
                self.sampleFits[-1].y.append(float(tokens[1]))
                self.sampleFits[-1].yp.append(float(tokens[2]))
                tokens=(tokens[3])[1:-1].split(',')
                self.sampleFits[-1].yl.append(float(tokens[0]))
                self.sampleFits[-1].yu.append(float(tokens[1]))

        fh.close()
