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

    unitDictionary = {\
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

    ROLE_TIME = 1010
    ROLE_MEASUREMENT = 1011
    ROLE_LABEL = 1012

    def __init__(self,tokens):
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
        typeString = ""
        if self.varType == PKPDVariable.TYPE_NUMERIC:
            typeString = "numeric"
        elif self.varType == PKPDVariable.TYPE_TEXT:
            typeString = "text"

        displayString=self.displayString.replace("%%","%%%%")

        roleString = ""
        if self.role == PKPDVariable.ROLE_TIME:
            roleString = "time"
        elif self.role == PKPDVariable.ROLE_MEASUREMENT:
            roleString = "measurement"
        elif self.role == PKPDVariable.ROLE_LABEL:
            roleString = "label"
        fh.write("%s ; %s ; %s[%s] ; %s ; %s\n"%(self.varName,self.units._toString(),
                                               typeString,displayString,roleString,self.comment))

class PKPDDose:
    TYPE_BOLUS = 1
    TYPE_INFUSION = 2

    def __init__(self,tokens):
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
        self.doseAmount = float(doseTokens[2].strip().lower().split("=")[1])

        # Get units
        unitString = tokens[2].strip()
        self.units = PKPDUnit(unitString)
        if not self.units.unit:
            raise Exception("Unrecognized unit: %s"%unitString)
        if not self.units.isWeight():
            raise Exception("After normalization, the dose must be a weight")

        # Get normalization
        self.normalization = tokens[3].strip()

    def _printToStream(self,fh):
        doseString = ""
        if self.doseType == PKPDDose.TYPE_BOLUS:
            doseString = "bolus t=%f"%self.t0
        elif self.doseType == PKPDDose.TYPE_INFUSION:
            doseString = "infusion t=%f...%f"%(self.t0,self.tF)
        fh.write("%s ; %s d=%f; %s ; %s\n"%\
                 (self.varName,doseString, self.doseAmount, self.units._toString(), self.normalization))

class PKPDSample:
    def parseLine(self,tokens,variableDict,doseDict):
        # FemaleRat1; dose=Dose1; weight=207

        # Keep a pointer to variableDict and doseDict
        self.variableDictPtr = variableDict
        self.doseDictPtr = doseDict

        # Get name
        self.varName = tokens[0].strip()

        # Get dose
        doseName = tokens[1].split('=')[1].strip()
        if doseName in doseDict:
            self.doseName = doseName
            self.dose = doseDict[doseName]
        else:
            raise Exception("Unrecognized dose %s"%doseName)

        # Get rest of variables
        self.descriptors = {}
        for n in range(2,len(tokens)):
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
                raise Exception("Unrocgnized variable %s"%varName)

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

    def _printToStream(self,fh):
        descriptorString = ""
        for key, value in self.descriptors.iteritems():
            descriptorString +="; %s=%s"%(key,value)
        fh.write("%s; dose=%s %s\n"%(self.varName,self.doseName,descriptorString))

    def _printMeasurements(self,fh):
        patternString = ""
        for n in range(0,len(self.measurementPattern)):
            patternString += "; %s"%self.measurementPattern[n]
        fh.write("%s %s\n"%(self.varName,patternString))
        if len(self.measurementPattern)>0:
            aux=getattr(self,"measurement_%s"%self.measurementPattern[0])
            N = len(aux)
        for i in range(0,N):
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

    def load(self, fnExperiment=""):
        if fnExperiment!="":
            self.fnPKPD.set(fnExperiment)
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
                self.variables[varname] = PKPDVariable(tokens)

            elif state==PKPDExperiment.READING_DOSES:
                tokens = line.split(';')
                if len(tokens)!=4:
                    print("Skipping dose: ",line)
                    continue
                dosename = tokens[0].strip()
                self.doses[dosename] = PKPDDose(tokens)

            elif state==PKPDExperiment.READING_SAMPLES:
                tokens = line.split(';')
                if len(tokens)<2:
                    print("Skipping sample: ",line)
                    continue
                samplename = tokens[0].strip()
                self.samples[samplename] = PKPDSample()
                self.samples[samplename].parseLine(tokens,self.variables, self.doses)

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

    def setExperiment(self, experiment):
        self.experiment = experiment
        self.fnExperiment = experiment.fnPKPD

    def setXVar(self, x):
        if not x in self.experiment.variables:
            raise Exception("Cannot find %s as a variable in the experiment"%x)
        self.xName = x
        self.xRange = self.experiment.getRange(x)

    def setX(self, xString):
        self.xString = xString

    def setY(self, yString):
        self.yString = yString

    def setYVar(self, y):
        if not y in self.experiment.variables:
            raise Exception("Cannot find %s as a variable in the experiment"%x)
        self.yName = y
        self.yRange = self.experiment.getRange(y)

    def xyToFloat(self):
        if len(self.xString)!=len(self.yString):
            raise Exception("Cannot model with different X and Y lengths")
        self.x = []
        self.y = []
        for n in range(0,len(self.xString)):
            if self.xString[n]!="NA" and self.yString[n]!="NA":
                self.x.append(float(self.xString[n]))
                self.y.append(float(self.yString[n]))

    def forwardModel(self):
        pass

    def getNumberOfParameters(self):
        pass

    def printSetup(self):
        pass

class PKPDExponentialModel(PKPDModel):
    def forwardModel(self, parameters):
        self.yPredicted = []
        for x in self.x:
            y = 0
            for k in range(0,self.Nexp):
                ck = parameters[2*k]
                lk = parameters[2*k+1]
                y += ck*math.exp(-lk*x)
            self.yPredicted.append(y)
        return self.yPredicted

    def setNexp(self, Nexp):
        self.Nexp = Nexp

    def getNumberOfParameters(self):
        return 2*self.Nexp

    def parseBounds(self,inputString,msg):
        if ';' in inputString:
            tokens=inputString.split(';')
            if len(tokens)!=self.Nexp:
                raise Exception("The number of intervals for %s does not match the number of exponential terms"%msg)
        else:
            tokens = [inputString]*self.Nexp
        return eval("["+",".join(tokens)+"]")

    def setBounds(self, cBounds, lambdaBounds):
        y0=self.yRange[1]
        yF=self.yRange[0]
        t0=self.xRange[0]
        tF=self.xRange[1]
        # y0 = c * exp(-lambda*t0)
        # yF = c * exp(-lambda*tF)
        meanLambda = math.log(yF/y0)/(t0-tF)
        meanC = y0*math.exp(meanLambda*t0)

        if lambdaBounds=="" or lambdaBounds is None:
            self.lambdaBounds = [(0.2*meanLambda,5*meanLambda)]*self.Nexp
        else:
            self.lambdaBounds=self.parseBounds(lambdaBounds,"lambda")

        if cBounds=="" or cBounds is None:
            self.cBounds = [(0.2*meanC,5*meanC)]*self.Nexp
        else:
            self.cBounds=self.parseBounds(cBounds,"amplitude")

    def getBounds(self):
        bounds = []
        for i in range(0,self.Nexp):
            bounds.append(self.cBounds[i])
            bounds.append(self.lambdaBounds[i])
        return bounds

    def printSetup(self):
        print("Model: y=sum_i c_i*exp(-lambda_i * x)")
        print("Number of exponentials: "+str(self.Nexp))
        print("Amplitude (c) bounds: "+str(self.cBounds))
        print("Lambda bounds: "+str(self.lambdaBounds))

class PKPDDEOptimizer:
    def __init__(self,model,fitType,goalFunction="RMSE"):
        self.model = model

        self.yTarget = np.array(model.y, dtype=np.float32)
        if fitType=="linear":
            self.takeYLogs = False
        elif fitType=="log":
            self.yTarget = np.log10(self.yTarget)
            self.takeYLogs = True

        if goalFunction=="RMSE":
            self.goalFunction = self.goalRMSE
        else:
            raise Exception("Unknown goal function")

    def goalRMSE(self,parameters):
        yPredicted = np.array(self.model.forwardModel(parameters),dtype=np.float32)
        if self.takeYLogs:
            yPredicted = np.log10(yPredicted)
        rmse = math.sqrt(((self.yTarget - yPredicted) ** 2).mean())
        print(rmse)
        return rmse

    def optimize(self):
        from scipy.optimize import differential_evolution
        self.optimum = differential_evolution(self.goalFunction, self.model.getBounds())
        return self.optimum

