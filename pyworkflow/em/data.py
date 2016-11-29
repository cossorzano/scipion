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

import copy
import json
import math
import numpy as np

from pyworkflow.em.pkpd_units import PKPDUnit, convertUnits
from pyworkflow.object import *
from pyworkflow.utils.path import writeMD5, verifyMD5


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
        self.units = None

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

    def isNumeric(self):
        return self.varType == self.TYPE_NUMERIC

    def isLabel(self):
        return self.role == self.ROLE_LABEL

    def isMeasurement(self):
        return self.role == self.ROLE_MEASUREMENT

    def isTime(self):
        return self.role == self.ROLE_TIME

    def getLabel(self):
        return "%s [%s]" % (self.varName, self.getUnitsString())


class PKPDDose:
    TYPE_BOLUS = 1
    TYPE_REPEATED_BOLUS = 2
    TYPE_INFUSION = 3

    def __init__(self):
        self.varName = None
        self.doseType = None
        self.doseAmount = None
        self.t0 = None
        self.tF = None
        self.every = None
        self.tunits = None
        self.dunits = None

    def parseTokens(self,tokens):
        # Dose1; bolus t=0 d=60*$(weight)/1000; min; mg
        # Dose1; repeated_bolus t=0:8:48 d=60*$(weight)/1000; h; mg
        # Dose1; infusion t=0:59 d=1; min; mg

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
        elif doseTypeString=="repeated_bolus":
            self.doseType = PKPDDose.TYPE_REPEATED_BOLUS
            timeTokens = timeString.split(":")
            self.t0 = float(timeTokens[0].strip())
            self.every = float(timeTokens[1].strip())
            self.tF = float(timeTokens[2].strip())
        elif doseTypeString=="infusion":
            self.doseType = PKPDDose.TYPE_INFUSION
            timeTokens = timeString.split(":")
            self.t0 = float(timeTokens[0].strip())
            self.tF = float(timeTokens[1].strip())
        else:
            raise Exception("Unrecognized dose type %s"%doseTypeString)
        self.doseAmount = doseTokens[2].strip().lower().split("=")[1]

        # Get time units
        unitString = tokens[2].strip()
        self.tunits = PKPDUnit(unitString)
        if not self.tunits.unit:
            raise Exception("Unrecognized unit: %s"%unitString)
        if not self.tunits.isTime():
            raise Exception("Time unit is not valid")

        # Get dose units
        unitString = tokens[3].strip()
        self.dunits = PKPDUnit(unitString)
        if not self.dunits.unit:
            raise Exception("Unrecognized unit: %s"%unitString)
        if not self.dunits.isWeight():
            raise Exception("After normalization, the dose must be a weight")

    def _printToStream(self,fh):
        fh.write("%s ; %s d=%s; %s; %s\n" % (self.varName,
                                              self.getDoseString(),
                                              self.doseAmount,
                                              self.tunits._toString(),
                                              self.dunits._toString()))

    def getDoseString(self):
        if self.doseType == PKPDDose.TYPE_BOLUS:
            doseString = "bolus t=%f" % self.t0
        elif self.doseType == PKPDDose.TYPE_REPEATED_BOLUS:
            doseString = "repeated_bolus t=%f:%f:%f" % (self.t0, self.every, self.tF)
        elif self.doseType == PKPDDose.TYPE_INFUSION:
            doseString = "infusion t=%f:%f" % (self.t0, self.tF)
        else:
            doseString = ""
        return doseString

    def changeTimeUnitsToMinutes(self):
        if self.doseType == PKPDDose.TYPE_BOLUS:
            self.t0 *= 60
        elif self.doseType == PKPDDose.TYPE_REPEATED_BOLUS:
            self.t0 *= 60
            self.tF *= 60
            self.every *= 60
        elif self.doseType == PKPDDose.TYPE_INFUSION:
            self.t0 *= 60
            self.tF *= 60

    def getDoseAt(self,t0,dt=0.5):
        """Dose between t0<=t<t0+dt, t0 is in minutes"""
        t1=t0+dt
        if self.doseType == PKPDDose.TYPE_BOLUS:
            if t0<=self.t0 and self.t0<t1:
                return self.doseAmount
            else:
                return 0.0
        elif self.doseType == PKPDDose.TYPE_REPEATED_BOLUS:
            doseAmount=0
            for t in np.arange(self.t0,self.tF,self.every):
                if t0<=t and t<t1:
                    doseAmount+=self.doseAmount
            return doseAmount
        elif self.doseType == PKPDDose.TYPE_INFUSION:
            if t0>self.tF or t1<self.t0:
                return 0.0
            else:
                tLeft=max(t0,self.t0)
                tRight=min(t1,self.tF)
                return self.doseAmount*(tRight-tLeft)

    def isDoseABolus(self):
        if self.doseType != PKPDDose.TYPE_BOLUS:
            return False
        if self.t0 != 0:
            return False
        return True

    def getTUnitsString(self):
        return self.tunits._toString()

    def getDUnitsString(self):
        return self.dunits._toString()


def createDeltaDose(doseAmount,t=0,dunits="mg"):
    dose = PKPDDose()
    dose.varName = "Bolus"
    dose.doseType = PKPDDose.TYPE_BOLUS
    dose.doseAmount = doseAmount
    dose.t0 = t
    dose.tunits = PKPDUnit("min")
    dose.dunits = PKPDUnit(dunits)
    return dose

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

    def getSampleName(self):
        return self.varName

    def interpretDose(self):
        self.parsedDoseList = []
        firstUnit = None
        for doseName in sorted(self.doseList):
            dose = copy.copy(self.doseDictPtr[doseName])
            dose.doseAmount = self.evaluateExpression(dose.doseAmount)
            dose.changeTimeUnitsToMinutes()
            if firstUnit==None:
                firstUnit = dose.dunits.unit
            else:
                dose.doseAmount = convertUnits(dose.doseAmount, dose.dunits.unit, firstUnit)
            self.parsedDoseList.append(dose)
        if len(self.parsedDoseList)==0:
            raise Exception("Cannot find any useful dose")

    def isDoseABolus(self):
        if len(self.parsedDoseList)!=1:
            return False
        return self.parsedDoseList[0].isDoseABolus()

    def getDoseAt(self,t0,dt=0.5):
        doseAmount = 0.0
        for dose in self.parsedDoseList:
            doseAmount += dose.getDoseAt(t0,dt)
        return doseAmount

    def getCumulatedDose(self,t0,tF):
        return self.getDoseAt(t0,tF-t0)

    def getDoseUnits(self):
        return self.parsedDoseList[0].dunits.unit

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
            if tokens[n]=="NA" or tokens[n]=="ULOQ" or tokens[n]=="LLOQ":
                ok = (self.variableDictPtr[varName].role != PKPDVariable.ROLE_TIME)
            if ok:
                exec("self.measurement_%s.append('%s')"%(varName,tokens[n]))
            else:
                raise Exception("Time measurements cannot be NA")

    def addMeasurementColumn(self,varName,values):
        self.measurementPattern.append(varName)
        setattr(self, "measurement_%s"%varName, [])
        for value in values:
            exec("self.measurement_%s.append('%s')"%(varName,str(value)))

    def getNumberOfVariables(self):
        return len(self.measurementPattern)

    def getNumberOfMeasurements(self):
        return len(getattr(self,"measurement_%s"%self.measurementPattern[0]))

    def _printToStream(self,fh):
        descriptorString = ""
        for key in sorted(self.descriptors.keys()):
            descriptorString +="; %s=%s"%(key,self.descriptors[key])
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
            aux = [x for x in aux if x != "NA" and x!="LLOQ" and x!="ULOQ"]
            x = np.asarray(aux, dtype=np.double)
            return [x.min(),x.max()]

    def getValues(self, varName):
        if type(varName)==list:
            retval=[]
            for vName in varName:
                if vName not in self.measurementPattern:
                    retval.append(None)
                else:
                    retval.append(getattr(self,"measurement_%s"%vName))
            return retval
        else:
            if varName not in self.measurementPattern:
                return None
            else:
                return getattr(self,"measurement_%s"%varName)

    def setValues(self, varName, varValues):
        setattr(self,"measurement_%s"%varName,varValues)

    def getXYValues(self,varNameX,varNameY):
        xl = []
        yl = []
        xs = self.getValues(varNameX)
        if type(varNameY)==list:
            ys = self.getValues(varNameY)
            for ysi in ys:
                xPartial =[]
                yPartial = []
                for x, y in izip(xs, ysi):
                    if x != "NA" and x!="LLOQ" and y!="ULOQ" and y != "NA" and y!= "LLOQ" and y!="ULOQ":
                        xPartial.append(float(x))
                        yPartial.append(float(y))
                xl.append(xPartial)
                yl.append(yPartial)
        else:
            ys = self.getValues(varNameY)
            for x, y in izip(xs, ys):
                if x != "NA" and x!="LLOQ" and y!="ULOQ" and y != "NA" and y!= "LLOQ" and y!="ULOQ":
                    xl.append(float(x))
                    yl.append(float(y))
        return xl, yl

    def getSampleMeasurements(self):
        return [PKPDSampleMeasurement(self,n) for n in range(0,self.getNumberOfMeasurements())]

    def substituteValuesInExpression(self, expression, prefix=""):
        expressionPython = copy.copy(expression)
        for key, variable in self.variableDictPtr.iteritems():
            if key in self.descriptors:
                value = self.descriptors[key]
                if value=="NA" or value=="LLOQ" or value=="ULOQ":
                    expressionPython="None"
                    break
                else:
                    if variable.varType == PKPDVariable.TYPE_NUMERIC:
                        expressionPython = expressionPython.replace("$%s(%s)"%(prefix,key),"%f"%float(value))
                    else:
                        expressionPython = expressionPython.replace("$%s(%s)"%(prefix,key),"'%s'"%value)
        return expressionPython

    def evaluateExpression(self, expression, prefix=""):
        expressionPython=self.substituteValuesInExpression(expression,prefix)
        return eval(expressionPython, {"__builtins__" : {"True": True, "False": False, "None": None} }, {})


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
        self.infoStr = String()
        self.general = {}
        self.variables = {}
        self.samples = {}
        self.doses = {}

    def __str__(self):
        if not self.infoStr.hasValue():
            self.load()
            self.infoStr.set("variables: %d, samples: %d"
                             % (len(self.variables), len(self.samples)))
        return self.infoStr.get()

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
                if len(tokens)!=4:
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
        self.infoStr.set("variables: %d, samples: %d" % (len(self.variables),
                                                         len(self.samples)))

    def _printToStream(self,fh):
        fh.write("[EXPERIMENT] ===========================\n")
        for key, value in self.general.iteritems():
            fh.write("%s = %s\n"%(key,value))
        fh.write("\n")

        fh.write("[VARIABLES] ============================\n")
        for key in sorted(self.variables.keys()):
            self.variables[key]._printToStream(fh)
        fh.write("\n")

        fh.write("[DOSES] ================================\n")
        for key in sorted(self.doses.keys()):
            self.doses[key]._printToStream(fh)
        fh.write("\n")

        fh.write("[SAMPLES] ================================\n")
        for key in sorted(self.samples.keys()):
            self.samples[key]._printToStream(fh)
        fh.write("\n")

        fh.write("[MEASUREMENTS] ===========================\n")
        for key in sorted(self.samples.keys()):
            self.samples[key]._printMeasurements(fh)
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

    def sampleSummary(self):
        summary=[]
        for varName, var in self.variables.iteritems():
            if var.role == PKPDVariable.ROLE_LABEL:
                toAdd = varName+": "
                if var.varType==PKPDVariable.TYPE_NUMERIC:
                    listOfValues=[]
                else:
                    listOfValues={}
                for sampleName, sample in self.samples.iteritems():
                    value = sample.descriptors[varName]
                    if var.varType==PKPDVariable.TYPE_NUMERIC:
                        listOfValues.append(float(value))
                    else:
                        if value in listOfValues:
                            listOfValues[value]+=1
                        else:
                            listOfValues[value]=1
                if var.varType==PKPDVariable.TYPE_NUMERIC:
                    listOfValuesNp = np.array(listOfValues)
                    toAdd += " mean=%f std=%f 5%%=%f 25%%=%f 50%%=%f 75%%=%f 95%%=%f"%\
                             (np.mean(listOfValuesNp),np.std(listOfValuesNp),np.percentile(listOfValuesNp,5),\
                              np.percentile(listOfValuesNp,25),np.percentile(listOfValuesNp,50),\
                              np.percentile(listOfValuesNp,75),np.percentile(listOfValuesNp,95))
                else:
                    for value in listOfValues:
                        toAdd += value + "(" + str(listOfValues[value]) + ") "
                summary.append(toAdd)
        return summary

    def getVarUnits(self,varName):
        if varName in self.variables:
            return self.variables[varName].units.unit
        else:
            return PKPDUnit.UNIT_NONE

    def addParameterToSample(self, sampleName, varName, varUnits, varDescr, varValue):
        if not varName in self.variables:
            varX = PKPDVariable()
            varX.varName = varName
            varX.varType = PKPDVariable.TYPE_NUMERIC
            varX.displayString = "%f"
            varX.role = PKPDVariable.ROLE_LABEL
            varX.comment = varDescr
            varX.units = PKPDUnit()
            varX.units.unit = varUnits
            self.variables[varName] = varX
        if sampleName in self.samples:
            sample = self.samples[sampleName]
            sample.descriptors[varName] = varValue

    def getSubGroup(self,condition):
        if condition=="":
            return self.samples
        samplesSubGroup = {}
        for sampleName, sample in self.samples.iteritems():
            if sample.evaluateExpression(condition):
                samplesSubGroup[sampleName] = sample
        return samplesSubGroup

    def getSubGroupLabels(self,condition,labelName):
        subgroupLabels = []
        for sampleName, sample in self.samples.iteritems():
            if condition!="" and sample.evaluateExpression(condition) or condition=="":
                subgroupLabels.append(sample.descriptors[labelName])
        return subgroupLabels

    def getNonBolusDoses(self):
        nonBolusList = []
        for sampleName, sample in self.samples.iteritems():
            sample.interpretDose()
            if not sample.isDoseABolus():
                nonBolusList.append(sampleName)
        return nonBolusList

class PKPDModelBase:
    def __init__(self):
        self.fnExperiment = None
        self.parameters = None
        self.parameterUnits = None
        self.xName = None
        self.yName = None

    def setExperiment(self, experiment):
        self.experiment = experiment
        if experiment!=None:
            self.fnExperiment = experiment.fnPKPD

    def setXVar(self, x):
        if not x in self.experiment.variables:
            raise Exception("Cannot find %s as a variable in the experiment"%x)
        self.xName = x
        self.xRange = self.experiment.getRange(x)

    def setYVar(self, y):
        if type(y)==list:
            self.yName = []
            self.yRange = []
            for yi in y:
                if not yi in self.experiment.variables:
                    raise Exception("Cannot find %s as a variable in the experiment"%yi)
                self.yName.append(yi)
                self.yRange.append(self.experiment.getRange(yi))
        else:
            if not y in self.experiment.variables:
                raise Exception("Cannot find %s as a variable in the experiment"%y)
            self.yName = y
            self.yRange = self.experiment.getRange(y)

    def setXYValues(self, x, y):
        self.x = np.array(x)
        if type(y[0])!=list and type(y[0])!=np.ndarray:
            self.y = np.array(y)
            idx = np.logical_and(np.isfinite(self.x), np.isfinite(self.y))
            self.x = self.x[idx]
            self.y = self.y[idx]
            self.ylog = [math.log10(yi) if yi>0 else float("inf") for yi in self.y]
        else:
            self.y = [np.array(yi) for yi in y]
            self.ylog = [np.log10(yi) for yi in self.y]

    def getNumberOfParameters(self):
        return len(self.getParameterNames())

    def getDescription(self):
        pass

    def getParameterNames(self):
        pass

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form %s'%self.getModelEquation()]*self.getNumberOfParameters()

    def calculateParameterUnits(self,sample):
        pass

    def setParameters(self, parameters):
        self.parameters = parameters

class PKPDModelBase2(PKPDModelBase):
    def __init__(self):
        PKPDModelBase.__init__(self)
        self.bounds = None

    def forwardModel(self, parameters, x=None):
        pass

    def printSetup(self):
        print("Model: %s"%self.getModelEquation())
        print("Variables: "+str(self.getParameterNames()))
        print("Bounds: "+str(self.getBounds()))

    def setSample(self, sample):
        self.sample = sample
        self.Dunits = sample.getDoseUnits()

    def getEquation(self):
        pass

    def getModelEquation(self):
        pass

    def getParameterDescriptions(self):
        return ['Automatically fitted model of the form %s'%self.getModelEquation()]*self.getNumberOfParameters()
        pass

    def areParametersSignificant(self, lowerBound, upperBound):
        """
        :param lowerBound and upperBound: a numpy array of parameters
        :return: a list of string with "True", "False", "NA", "Suspicious"
        """
        pass

    def areParametersValid(self, p):
        pass

    def setBounds(self, boundsString):
        self.bounds = None
        if boundsString!="" and boundsString!=None:
            tokens=boundsString.split(';')
            if len(tokens)!=self.getNumberOfParameters():
                raise Exception("The number of bound intervals does not match the number of exponential terms")
            self.bounds=[]
            for token in tokens:
                values = token.strip().split(',')
                self.bounds.append((float(values[0][1:]),float(values[1][:-1])))

    def getBounds(self):
        return self.bounds

    def setConfidenceInterval(self,lowerBound,upperBound):
        yPredictedBackup = copy.copy(self.yPredicted)
        self.yPredictedUpper = copy.copy(self.yPredicted)
        self.yPredictedLower = copy.copy(self.yPredicted)
        for i in range(0,int(math.pow(2,self.getNumberOfParameters()))):
            pattern = ("{0:0%db}"%(self.getNumberOfParameters())).format(i)
            p = np.where(np.array(list(pattern))=="1",upperBound,lowerBound)
            p = p.astype(np.float)
            if not self.areParametersValid(p):
                continue
            y = self.forwardModel(p)
            if type(y[0])!=list and type(y[0])!=np.ndarray:
                for n in range(len(y)):
                    if y[n]<self.yPredictedLower[n]:
                        if y[n]<0:
                            self.yPredictedLower[n]=0
                        else:
                            self.yPredictedLower[n]=y[n]
                    if y[n]>self.yPredictedUpper[n]:
                        self.yPredictedUpper[n]=y[n]
            else:
                for j in range(len(y)):
                    yj=y[j]
                    for n in range(len(yj)):
                        if yj[n]<self.yPredictedLower[j][n]:
                            if yj[n]<0:
                                self.yPredictedLower[j][n]=0
                            else:
                                self.yPredictedLower[j][n]=yj[n]
                        if yj[n]>self.yPredictedUpper[j][n]:
                            self.yPredictedUpper[j][n]=yj[n]
        self.yPredicted = yPredictedBackup

    def setConfidenceIntervalNA(self):
        if type(self.yPredicted[0])!=np.ndarray:
            self.yPredictedUpper = ["NA"]*len(self.yPredicted)
            self.yPredictedLower = ["NA"]*len(self.yPredicted)
        else:
            self.yPredictedUpper = []
            self.yPredictedLower = []
            for y in self.yPredicted:
                self.yPredictedUpper.append(["NA"]*len(y))
                self.yPredictedLower.append(["NA"]*len(y))

class PKPDModel(PKPDModelBase2):
    def prepare(self):
        pass

class PKPDODEModel(PKPDModelBase2):
    def __init__(self):
        PKPDModelBase2.__init__(self)
        self.t0 = None # (min)
        self.tF = None # (min)
        self.deltaT = 0.25 # (min)
        self.drugSource = None
        # self.show = False

    def F(self, t, y):
        return 0

    def G(self, t, dD):
        return 0

    def imposeConstraints(self, yt):
        pass

    def getResponseDimension(self):
        return None

    def getStateDimension(self):
        return None

    def forwardModel(self, parameters, x=None):
        self.parameters = parameters

        # Simulate the system response
        t = self.t0
        Nsamples = int(math.ceil((self.tF-self.t0)/self.deltaT))+1
        if self.getStateDimension()>1:
            yt = np.zeros(self.getStateDimension(),np.double)
            Yt = np.zeros((Nsamples,self.getStateDimension()),np.double)
        else:
            yt = 0.0
            Yt = np.zeros(Nsamples)
        Xt = np.zeros(Yt.shape[0])
        delta_2 = 0.5*self.deltaT
        K = self.deltaT/3
        for i in range(0,Nsamples):
            t = self.t0 + i*self.deltaT # More accurate than t+= self.deltaT
            Xt[i]=t

            # Internal evolution
            # Runge Kutta's 4th order (http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html)
            k1 = self.F(t,yt)
            k2 = self.F(t+delta_2,yt+k1*delta_2)
            k3 = self.F(t+delta_2,yt+k2*delta_2)
            k4 = self.F(t+self.deltaT,yt+k3*self.deltaT)

            # External driving: Drug component
            dD = self.drugSource.getAmountReleasedBetween(t,t+self.deltaT)
            dyD = self.G(t, dD)

            # Update state
            yt += (0.5*(k1+k4)+k2+k3)*K+dyD

            # Make sure it makes sense
            self.imposeConstraints(yt)

            # if self.show:
            #     print("t=%f dD=%s dyD=%s dy=%s"%(t,str(dD),str(dyD),str((0.5*(k1+k4)+k2+k3)*K)))

            # Keep this result and go to next iteration
            if self.getStateDimension()>1:
                Yt[i,:]=yt
            else:
                Yt[i]=yt

        # Get the values at x
        if x==None:
            x = self.x

        if type(x[0])!=list and type(x[0])!=np.ndarray:
            if self.getStateDimension()==1:
                self.yPredicted = np.interp(x,Xt,Yt)
            else:
                self.yPredicted = []
                for j in range(0,self.getResponseDimension()):
                    self.yPredicted.append(np.interp(x,Xt,Yt[:,j]))
        else:
            self.yPredicted = []
            for j in range(0,self.getResponseDimension()):
                self.yPredicted.append(np.interp(x[j],Xt,Yt[:,j]))
        return self.yPredicted

class PKPDOptimizer:
    def __init__(self,model,fitType,goalFunction="RMSE"):
        self.model = model

        if type(model.y[0])!=list and type(model.y[0])!=np.ndarray:
            self.yTarget = np.array(model.y, dtype=np.float32)
            self.yTargetLogs = [np.array([math.log10(yi) if np.isfinite(yi) and yi>0 else float("inf") for yi in self.yTarget],dtype=np.float32)]
        else:
            self.yTarget = [np.array(yi, dtype=np.float32) for yi in model.y]
            self.yTargetLogs = [np.log10(yi) for yi in self.yTarget]

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

        self.bounds = model.getBounds()

        if goalFunction=="RMSE":
            self.goalFunction = self.goalRMSE
        else:
            raise Exception("Unknown goal function")

        self.verbose = 1

    def inBounds(self,parameters):
        if self.bounds==None or len(self.bounds)!=len(parameters):
            return True
        for n in range(0,len(parameters)):
            if parameters[n]<self.bounds[n][0] or parameters[n]>self.bounds[n][1]:
                return False
        return True

    def getResiduals(self,parameters):
        if not self.inBounds(parameters):
            if type(self.yTarget[0])==list or type(self.yTarget[0])==np.ndarray:
                allDiffs = None
                for yTarget in self.yTarget:
                    diff = 1e38*np.ones(yTarget.shape)
                    if allDiffs==None:
                        allDiffs = diff
                    else:
                        allDiffs = np.concatenate([allDiffs, diff])
                return allDiffs
            else:
                return 1e38*np.ones(self.yTarget.shape)
        yPredicted = self.model.forwardModel(parameters)
        if type(yPredicted[0])!=list and type(yPredicted[0])!=np.ndarray:
            yPredicted = [np.array(yPredicted,dtype=np.float32)]
        else:
            yPredicted = [np.array(yi,dtype=np.float32) for yi in yPredicted]

        allDiffs = None
        for y, yTarget in izip(yPredicted,self.yTarget):
            if self.takeYLogs:
                idx = np.array([np.isfinite(yi) and yi>=1e-20 for yi in y])
                nonIdx = np.logical_not(idx)
                y[idx] = np.log10(y[idx])
                y[nonIdx] = -100
            diff = yTarget - y
            if self.takeRelative:
                diff = diff/yTarget
            if allDiffs==None:
                allDiffs = diff
            else:
                allDiffs = np.concatenate([allDiffs, diff])
            idx = np.logical_not(np.isfinite(allDiffs))
            allDiffs[idx]=1e38
        return allDiffs[np.isfinite(allDiffs)]

    def goalRMSE(self,parameters):
        e = self.getResiduals(parameters)
        rmse = math.sqrt(np.power(e,2).mean())
        return rmse

    def _evaluateQuality(self, x, y, yp):
        # Spiess and Neumeyer, BMC Pharmacology 2010, 10:6
        if type(y[0])==list or type(y[0])==np.ndarray:
            self.e = None
            yToUse = None
            for yi, ypi in izip(y,yp):
                diff = yi - ypi
                if self.e==None:
                    self.e = diff
                    yToUse = np.asarray(yi)
                else:
                    self.e = np.concatenate([self.e, diff])
                    yToUse = np.concatenate([yToUse, yi])
        else:
            yToUse = np.asarray(y).flatten()
            self.e = yToUse-np.asarray(yp).flatten()

        self.R2 = (1-np.var(self.e)/np.var(yToUse))
        n=len(self.e) # Number of samples
        p=self.model.getNumberOfParameters()
        if n-p>0:
            self.R2adj = 1-self.R2*(n-1)/(n-p)*(1-self.R2)
        else:
            self.R2adj = -1
        logL = 0.5*(-n*(math.log(2*math.pi)+1-math.log(n)+math.log(np.sum(np.multiply(self.e,self.e)))))
        self.AIC = 2*p-2*logL
        if n-p-1>0:
            self.AICc = self.AIC+2*p*(p+1)/(n-p-1)
        else:
            self.AICc = -1
        self.BIC = p*math.log(n)-2*logL

    def _printFitting(self, x, y, yp):
        if type(y)==np.ndarray and y.size==0:
            return
        if type(y[0])!=list and type(y[0])!=np.ndarray and (type(yp[0])==list or type(yp[0])==np.ndarray):
            yp=yp[0]
        self._evaluateQuality(x, y, yp)
        if not hasattr(self.model,"model") or self.model.model.getResponseDimension()==1:
            for n in range(0,x.shape[0]):
                print("%f %f %f %f"%(x[n],y[n],yp[n],y[n]-yp[n]))
        else:
            for j in range(len(x)):
                print("Series %d ---------"%j)
                xj=np.asarray(x[j])
                yj=y[j]
                ypj=yp[j]
                for n in range(0,xj.shape[0]):
                    print("%f %f %f %f"%(xj[n],yj[n],ypj[n],yj[n]-ypj[n]))
        print("------------------------")
        print("Mean error = %f"%np.mean(self.e))
        print("Std error = %f"%np.std(self.e))
        print("R2 = %f"%self.R2)
        print("R2adj = %f"%self.R2adj)
        print("AIC = %f"%self.AIC)
        print("AICc(Recommended) = %f"%self.AICc)
        print("BIC = %f"%self.BIC)
        print("------------------------")

    def evaluateQuality(self):
        yPredicted = np.array(self.model.forwardModel(self.model.parameters),dtype=np.float32)
        self._evaluateQuality(self.model.x, self.model.y, yPredicted)

    def printFitting(self):
        yPredicted = self.model.forwardModel(self.model.parameters)
        print("==========================================")
        print("X     Y    Ypredicted  Error=Y-Ypredicted ")
        print("==========================================")
        self._printFitting(self.model.x, self.model.y, yPredicted)
        print("==================================================================")
        print("X    log10(Y)  log10(Ypredicted)  Error=log10(Y)-log10(Ypredicted)")
        print("==================================================================")
        if type(self.model.y[0])!=list and type(self.model.y[0])!=np.ndarray:
            ylog = np.array([math.log10(yi) if np.isfinite(yi) and yi>0 else float("inf") for yi in self.model.y])
            yplog = np.array([math.log10(ypi) if np.isfinite(ypi) and ypi>0 else float("inf") for ypi in yPredicted])
            idx = np.array(np.where(np.logical_and(np.isfinite(ylog),np.isfinite(yplog))))
            self._printFitting(self.model.x[idx].ravel(), ylog[idx].ravel(), yplog[idx].ravel())
        else:
            logY = [np.log10(np.asarray(y)) for y in self.model.y]
            logYp = [np.log10(np.asarray(y)) for y in yPredicted]
            self._printFitting(self.model.x, logY, logYp)

class PKPDDEOptimizer(PKPDOptimizer):
    def optimize(self):
        from scipy.optimize import differential_evolution
        if self.verbose>0:
            print("Optimizing with Differential Evolution (DE), a global optimizer")
        self.optimum = differential_evolution(self.goalFunction, self.model.getBounds())
        if self.verbose>0:
            print("Best DE function value: "+str(self.optimum.fun))
            print("Best DE parameters: "+str(self.optimum.x))
        self.model.setParameters(self.optimum.x)
        if self.verbose>0:
            print(self.model.getEquation())
            self.printFitting()
            print(" ")
        return self.optimum

class PKPDLSOptimizer(PKPDOptimizer):
    def optimize(self):
        from scipy.optimize import leastsq
        if self.verbose>0:
            print("Optimizing with Least Squares (LS), a local optimizer")
            print("Initial parameters: "+str(self.model.parameters))
        self.optimum, self.cov_x, self.info, mesg, _ = leastsq(self.getResiduals, self.model.parameters, full_output=True)
        if self.verbose>0:
            print("Best LS function value: "+str(self.goalFunction(self.optimum)))
            print("Best LS parameters: "+str(self.optimum))
            print("Covariance matrix:")
            if self.cov_x!=None:
                self.cov_x *= np.var(self.info["fvec"])
                print(np.array_str(self.cov_x,max_line_width=120))
            else:
                print("Singular covariance matrix, at least one of the variables seems to be irrelevant")
        self.model.setParameters(self.optimum)
        if self.verbose>0:
            print(self.model.getEquation())
            self.printFitting()
            print(" ")
        return self.optimum

    def setConfidenceInterval(self,confidenceInterval):
        if self.cov_x!=None:
            from scipy.stats import norm
            nstd = norm.ppf(1-(1-confidenceInterval/100)/2)
            perr = np.sqrt(np.diag(self.cov_x))
            self.lowerBound = self.optimum-nstd*perr
            self.upperBound = self.optimum+nstd*perr

            self.significance = self.model.areParametersSignificant(self.lowerBound,self.upperBound)
            parameterNames = self.model.getParameterNames()
            print("Confidence intervals %f%% --------------------------"%confidenceInterval)
            print("ParameterName ParameterValue   ParameterConfidenceInterval  IsStatisticallySignificant")
            for n in range(0,len(self.optimum)):
                print("%s %f [%f,%f] %s"%(parameterNames[n],self.optimum[n],self.lowerBound[n],self.upperBound[n],self.significance[n]))

            self.model.setConfidenceInterval(self.lowerBound,self.upperBound)
        else:
            self.lowerBound=["NA"]*len(self.optimum)
            self.upperBound=["NA"]*len(self.optimum)
            self.significance=["NA"]*len(self.optimum)
            self.model.setConfidenceIntervalNA()

class PKPDSampleFit:
    READING_SAMPLEFITTINGS_NAME = 0
    READING_SAMPLEFITTINGS_MODELEQ = 1
    READING_SAMPLEFITTINGS_R2 = 2
    READING_SAMPLEFITTINGS_R2ADJ = 3
    READING_SAMPLEFITTINGS_AIC = 4
    READING_SAMPLEFITTINGS_AICc = 5
    READING_SAMPLEFITTINGS_BIC = 6
    READING_SAMPLEFITTINGS_PARAMETER_BOUNDS = 7
    READING_SAMPLEFITTINGS_SAMPLE_VALUES = 8

    def __init__(self):
        self.sampleName = ""

        # Lists with the sample fits values
        self.x = None
        self.y = None
        self.yp = None
        self.yl = None
        self.yu = None

        self.modelEquation = ""
        self.R2 = 0
        self.R2adj = 0
        self.AIC = 0
        self.AICc = 0
        self.BIC = 0
        self.parameters = None
        self.lowerBound = None
        self.upperBound = None
        self.significance = None

        self.multiOutputSeries = False

    def printForPopulation(self,fh,observations):
        outputStr = ""
        for parameter in self.parameters:
            outputStr += "%f "%parameter
        observations = np.vstack([observations, self.parameters])
        outputStr += " # %f %f %f %f %f"%(self.R2,self.R2adj,self.AIC,self.AICc,self.BIC)
        fh.write(outputStr+"\n")
        return observations

    def getBasicInfo(self):
        """ Return a string with some basic information of the fitting. """
        info = ""
        info += "Sample name: %s\n" % self.sampleName
        info += "Model: %s\n" % self.modelEquation
        info += "R2: %f\n" % self.R2
        info += "R2adj: %f\n" % self.R2adj
        info += "AIC: %f\n" % self.AIC
        info += "AICc(Recommended): %f\n" % self.AICc
        info += "BIC: %f\n" % self.BIC

        return info

    def _printToStream(self, fh):
        if self.significance == None:
            self.significance = ["Undetermined"]*len(self.parameters)
        fh.write(self.getBasicInfo())
        fh.write("Parameter lowerBound upperBound IsStatisticallySignificant -------\n")
        for parameter, lower, upper, significance in izip(self.parameters,self.lowerBound,self.upperBound,\
                                                          self.significance):
            fh.write("%f [%s,%s] %s\n"%(parameter,str(lower),str(upper),significance))
        fh.write("X   Y   Ypredicted [Ylower,Yupper] -------\n")
        if  any(isinstance(el, list) for el in self.x):
            for j in range(len(self.x)):
                fh.write("Series %d -----\n"%j)
                xj = self.x[j]
                yj = self.y[j]
                ypj = self.yp[j]
                ylj = self.yl[j]
                yuj = self.yu[j]
                for x,y,yp,yl,yu in izip(xj,yj,ypj,ylj,yuj):
                    fh.write("%f %s %s [%s,%s]\n"%(x,str(y),str(yp),str(yl),str(yu)))
        else:
            for x,y,yp,yl,yu in izip(self.x,self.y,self.yp,self.yl,self.yu):
                fh.write("%f %s %s [%s,%s]\n"%(x,str(y),str(yp),str(yl),str(yu)))
        fh.write("\n")

    def restartReadingState(self):
        self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_NAME

    def readFromLine(self, line):
        if self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_NAME:
            tokens = line.split(':')
            self.sampleName = tokens[1].strip()
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_MODELEQ

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_MODELEQ:
            tokens = line.split(':')
            self.modelEquation = tokens[1].strip()
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_R2

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_R2:
            tokens = line.split(':')
            self.R2 = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_R2ADJ

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_R2ADJ:
            tokens = line.split(':')
            self.R2adj = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_AIC

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_AIC:
            tokens = line.split(':')
            self.AIC = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_AICc

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_AICc:
            tokens = line.split(':')
            self.AICc = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_BIC

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_BIC:
            tokens = line.split(':')
            self.BIC = float(tokens[1])
            self.state = PKPDSampleFit.READING_SAMPLEFITTINGS_PARAMETER_BOUNDS

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_PARAMETER_BOUNDS:
            if line.startswith("Parameter lowerBound upperBound"):
                self.parameters=[]
                self.lowerBound=[]
                self.upperBound=[]
                self.significance=[]
            elif line.startswith("X   Y   Ypredicted"):
                self.state=PKPDSampleFit.READING_SAMPLEFITTINGS_SAMPLE_VALUES
                self.x=[]
                self.y=[]
                self.yp=[]
                self.yl=[]
                self.yu=[]
            else:
                tokens=line.split()
                self.parameters.append(float(tokens[0]))
                self.significance.append(tokens[2])
                tokens=(tokens[1])[1:-1].split(',')
                self.lowerBound.append(tokens[0])
                self.upperBound.append(tokens[1])

        elif self.state==PKPDSampleFit.READING_SAMPLEFITTINGS_SAMPLE_VALUES:
            tokens=line.split()
            if tokens[0].strip()=="Series":
                self.x.append([])
                self.y.append([])
                self.yp.append([])
                self.yl.append([])
                self.yu.append([])
                self.multiOutputSeries = True
            else:
                if self.multiOutputSeries:
                    self.x[-1].append(float(tokens[0]))
                    self.y[-1].append(float(tokens[1]))
                    self.yp[-1].append(float(tokens[2]))
                    tokens=(tokens[3])[1:-1].split(',')
                    self.yl[-1].append(tokens[0])
                    self.yu[-1].append(tokens[1])
                else:
                    self.x.append(float(tokens[0]))
                    self.y.append(float(tokens[1]))
                    self.yp.append(float(tokens[2]))
                    tokens=(tokens[3])[1:-1].split(',')
                    self.yl.append(tokens[0])
                    self.yu.append(tokens[1])

    def copyFromOptimizer(self,optimizer):
        self.R2 = optimizer.R2
        self.R2adj = optimizer.R2adj
        self.AIC = optimizer.AIC
        self.AICc = optimizer.AICc
        self.BIC = optimizer.BIC
        self.significance = optimizer.significance
        self.lowerBound = optimizer.lowerBound
        self.upperBound = optimizer.upperBound


class PKPDSampleFitBootstrap:
    READING_SAMPLEFITTINGS_NAME = 0
    READING_SAMPLEFITTINGS_XB = 1
    READING_SAMPLEFITTINGS_YB = 2
    READING_SAMPLEFITTINGS_PARAMETERS = 3

    def __init__(self):
        self.sampleName = ""
        self.R2 = []
        self.R2adj = []
        self.AIC = []
        self.AICc = []
        self.BIC = []
        self.parameters = None
        self.xB = []
        self.yB = []

    def _printSample(self,fh,n, n0=0):
        outputStr = "%d: "%(n+n0)
        for parameter in self.parameters[n,:]:
            outputStr += "%f "%parameter
        outputStr += " # %f %f %f %f %f"%(self.R2[n],self.R2adj[n],self.AIC[n],self.AICc[n],self.BIC[n])
        fh.write(outputStr+"\n")

    def printForPopulation(self,fh,observations):
        for n in range(0,self.parameters.shape[0]):
            self._printSample(fh,n,observations.shape[0])
        observations = np.vstack([observations, self.parameters])
        return observations

    def _printToStream(self,fh):
        fh.write("Sample name: %s\n"%self.sampleName)
        for n in range(0,self.parameters.shape[0]):
            fh.write("xB: %s\n"%self.xB[n])
            fh.write("yB: %s\n"%self.yB[n])

            outputStr = ""
            for parameter in self.parameters[n,:]:
                outputStr += "%f "%parameter
            outputStr += " # %f %f %f %f %f"%(self.R2[n],self.R2adj[n],self.AIC[n],self.AICc[n],self.BIC[n])
            fh.write(outputStr+"\n")
        fh.write("\n")

    def restartReadingState(self):
        self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_NAME

    def readFromLine(self, line):
        if self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_NAME:
            tokens = line.split(':')
            self.sampleName = tokens[1].strip()
            self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_XB

        elif self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_XB:
            tokens = line.split(':')
            self.xB.append(tokens[1].strip())
            self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_YB

        elif self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_YB:
            tokens = line.split(':')
            self.yB.append(tokens[1].strip())
            self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_PARAMETERS

        elif self.state==PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_PARAMETERS:
            tokens = line.split('#')
            tokensParameters = tokens[0].strip().split(' ')
            tokensQuality = tokens[1].strip().split(' ')
            if self.parameters == None:
                self.parameters = np.empty((0,len(tokensParameters)),np.double)
            self.parameters = np.vstack([self.parameters, [float(prm) for prm in tokensParameters]])

            self.R2.append(float(tokensQuality[0]))
            self.R2adj.append(float(tokensQuality[1]))
            self.AIC.append(float(tokensQuality[2]))
            self.AICc.append(float(tokensQuality[3]))
            self.BIC.append(float(tokensQuality[4]))

            self.state = PKPDSampleFitBootstrap.READING_SAMPLEFITTINGS_XB

    def copyFromOptimizer(self,optimizer):
        self.R2.append(optimizer.R2)
        self.R2adj.append(optimizer.R2adj)
        self.AIC.append(optimizer.AIC)
        self.AICc.append(optimizer.AICc)
        self.BIC.append(optimizer.BIC)


class PKPDFitting(EMObject):
    READING_FITTING_EXPERIMENT = 1
    READING_FITTING_PREDICTOR = 2
    READING_FITTING_PREDICTED = 3
    READING_FITTING_PREDICTED_LIST = 4
    READING_FITTING_MODEL = 5
    READING_POPULATION_HEADER = 6
    READING_POPULATION = 7
    READING_SAMPLEFITTINGS_BEGIN = 8
    READING_SAMPLEFITTINGS_CONTINUE = 9

    def __init__(self, cls="", **args):
        EMObject.__init__(self, **args)
        self.fnFitting = String()
        self.fnExperiment = String()
        self.predictor = None
        self.predicted = None
        self.modelDescription = ""
        self.modelParameters = []
        self.modelParameterUnits = []
        self.sampleFits = []
        self.summaryLines = []
        if cls=="":
            self.sampleFittingClass = "PKPDSampleFit"
        else:
            self.sampleFittingClass = cls

    def isPopulation(self):
        return self.fnFitting.get().endswith("bootstrapPopulation.pkpd")

    def write(self, fnFitting):
        fh=open(fnFitting,'w')
        self._printToStream(fh)
        fh.close()
        self.fnFitting.set(fnFitting)
        writeMD5(fnFitting)

    def getAllParameters(self):
        allParameters = np.empty((0,len(self.modelParameters)),np.double)
        for sampleFitting in self.sampleFits:
            allParameters = np.vstack([allParameters, sampleFitting.parameters])
        return allParameters

    def getStats(self, observations=None):
        if observations is None:
            observations = self.getAllParameters()
        mu=np.mean(observations,axis=0)
        C=np.cov(np.transpose(observations))
        sigma = np.sqrt(np.diag(C))
        R=np.corrcoef(np.transpose(observations))
        percentiles = np.percentile(observations,[0, 2.5, 25, 50, 75, 97.5, 100],axis=0)
        return mu, sigma, R, percentiles

    def _printToStream(self,fh):
        fh.write("[FITTING] ===========================\n")
        fh.write("Experiment: %s\n"%self.fnExperiment.get())
        fh.write("Predictor (X): ")
        self.predictor._printToStream(fh)
        fh.write("Predicted (Y): ")
        if type(self.predicted)==list:
            fh.write("Predicted list=%d\n"%len(self.predicted))
            for y in self.predicted:
                y._printToStream(fh)
        else:
            self.predicted._printToStream(fh)
        fh.write("Model: %s\n"%self.modelDescription)
        fh.write("\n")

        fh.write("[POPULATION PARAMETERS] =============\n")
        auxUnit = PKPDUnit()
        for paramName, paramUnits in izip(self.modelParameters, self.modelParameterUnits):
            auxUnit.unit = paramUnits
            fh.write("%s [%s] "%(paramName,auxUnit._toString()))
        fh.write(" # R2 R2adj AIC AICc BIC\n")
        observations = np.empty((0,len(self.modelParameters)),np.double)
        for sampleFitting in self.sampleFits:
            observations = sampleFitting.printForPopulation(fh,observations)
        fh.write("\n")

        mu=np.mean(observations,axis=0)
        if observations.shape[0]>2:
            C=np.cov(np.transpose(observations))
            R=np.corrcoef(np.transpose(observations))
            sigma = np.sqrt(np.diag(C))
            fh.write("Mean   parameters  = %s\n"%np.array_str(mu))
            fh.write("Median parameters  = %s\n"%np.array_str(np.median(observations,axis=0)))
            limits = np.percentile(observations,[2.5,97.5],axis=0)
            fh.write("Lower bound (2.5%%) = %s\n"%np.array_str(limits[0]))
            fh.write("Upper bound (97.5%%) = %s\n"%np.array_str(limits[1]))
            fh.write("Covariance matrix  =\n%s\n"%np.array_str(C,max_line_width=120))
            fh.write("Correlation matrix  =\n%s\n"%np.array_str(R,max_line_width=120))
        fh.write("\n")

        fh.write("[SAMPLE FITTINGS] ===================\n")
        for sampleFitting in self.sampleFits:
            sampleFitting._printToStream(fh)

    def load(self, fnFitting=None):
        fnFitting = str(fnFitting or self.fnFitting)
        fh = open(fnFitting)
        if not fh:
            raise Exception("Cannot open %s" % fnFitting)
        if not verifyMD5(fnFitting):
            raise Exception("The file %s has been modified since its creation" % fnFitting)
        self.fnFitting.set(fnFitting)

        auxUnit = PKPDUnit()
        for line in fh.readlines():
            line=line.strip()
            if line=="":
                if state==PKPDFitting.READING_SAMPLEFITTINGS_CONTINUE:
                     state=PKPDFitting.READING_SAMPLEFITTINGS_BEGIN
                continue
            if line.startswith('[') and line.endswith('='):
                section = line.split('=')[0].strip().lower()
                if section=="[fitting]":
                    state=PKPDFitting.READING_FITTING_EXPERIMENT
                    self.summaryLines.append(line)
                elif section=="[population parameters]":
                    state=PKPDFitting.READING_POPULATION_HEADER
                    self.summaryLines.append(line)
                elif section=="[sample fittings]":
                    state=PKPDFitting.READING_SAMPLEFITTINGS_BEGIN
                else:
                    print("Skipping: ",line)

            elif state==PKPDFitting.READING_FITTING_EXPERIMENT:
                tokens = line.split(':')
                self.fnExperiment.set(tokens[1].strip())
                state = PKPDFitting.READING_FITTING_PREDICTOR
                self.summaryLines.append(line)

            elif state==PKPDFitting.READING_FITTING_PREDICTOR:
                tokens = line.split(':')
                self.predictor = PKPDVariable()
                self.predictor.parseTokens(tokens[1].split(';'))
                state = PKPDFitting.READING_FITTING_PREDICTED
                self.summaryLines.append(line)

            elif state==PKPDFitting.READING_FITTING_PREDICTED:
                tokens = line.split(':')
                if (tokens[1].strip().startswith("Predicted list")):
                    state = PKPDFitting.READING_FITTING_PREDICTED_LIST
                    self.remainingPredicted = int(tokens[1].split('=')[1])
                    self.predicted = []
                    self.summaryLines.append(line)
                else:
                    self.predicted = PKPDVariable()
                    self.predicted.parseTokens(tokens[1].split(';'))
                    state = PKPDFitting.READING_FITTING_MODEL
                    self.summaryLines.append(line)

            elif state==PKPDFitting.READING_FITTING_PREDICTED_LIST:
                    self.summaryLines.append(line)
                    newVar = PKPDVariable()
                    newVar.parseTokens(line.split(';'))
                    self.predicted.append(newVar)
                    self.remainingPredicted -= 1
                    if self.remainingPredicted == 0:
                        state = PKPDFitting.READING_FITTING_MODEL

            elif state==PKPDFitting.READING_FITTING_MODEL:
                tokens = line.split(':')
                self.modelDescription = tokens[1].strip()
                state = PKPDFitting.READING_POPULATION
                self.summaryLines.append(line)
                self.summaryLines.append("\n")

            elif state==PKPDFitting.READING_POPULATION_HEADER:
                lineParts = line.split('#')
                tokens = lineParts[0].strip().split(' ')
                for i in range(0,len(tokens),2):
                    self.modelParameters.append(tokens[i])
                    self.modelParameterUnits.append(auxUnit._fromString(tokens[i+1][1:-1]))
                state = PKPDFitting.READING_POPULATION

            elif state==PKPDFitting.READING_POPULATION:
                self.summaryLines.append(line)

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_BEGIN:
                newSampleFit = eval("%s()"%self.sampleFittingClass)
                self.sampleFits.append(newSampleFit)
                self.sampleFits[-1].restartReadingState()
                self.sampleFits[-1].readFromLine(line)
                state = PKPDFitting.READING_SAMPLEFITTINGS_CONTINUE

            elif state==PKPDFitting.READING_SAMPLEFITTINGS_CONTINUE:
                self.sampleFits[-1].readFromLine(line)

        fh.close()

    def getSampleFit(self, sampleName):
        for sampleFit in self.sampleFits:
            if sampleFit.sampleName == sampleName:
                return sampleFit
        return None

    def loadExperiment(self):
        experiment = PKPDExperiment()
        experiment.load(self.fnExperiment.get())
        return experiment


class PKPDSampleSignalAnalysis:
    def __init__(self):
        self.sampleName = ""
        self.x = None
        self.y = None
        self.analysisVariables = []
        self.parameters = None
        self.lowerBound = None
        self.upperBound = None
        self.significance = None

    def _printToStream(self,fh):
        if self.significance == None:
            self.significance = ["Undetermined"]*len(self.parameters)
        fh.write("Sample name: %s\n"%self.sampleName)
        for varName, parameter in izip(self.analysisVariables,self.parameters):
            fh.write("%s = %f\n"%(varName,parameter))
        fh.write("\n")


class PKPDSignalAnalysis(EMObject):
    READING_EXPERIMENT = 1
    READING_PREDICTOR = 2
    READING_PREDICTED = 3
    READING_MODEL = 4
    READING_ANALYSIS_VARIABLES = 5
    READING_POPULATION_HEADER = 6
    READING_POPULATION = 7
    READING_SAMPLEANALYSIS_NAME = 8
    READING_SAMPLEANALYSIS_PARAMETERS = 9

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.fnAnalysis = String()
        self.fnExperiment = String()
        self.predictor = None
        self.predicted = None
        self.analysisDescription = ""
        self.analysisParameters = []
        self.sampleAnalyses = []
        self.summaryLines = []

    def write(self, fnAnalysis):
        fh=open(fnAnalysis,'w')
        self._printToStream(fh)
        fh.close()
        self.fnAnalysis.set(fnAnalysis)
        writeMD5(fnAnalysis)

    def _printToStream(self,fh):
        fh.write("[SIGNAL ANALYSIS] ===========================\n")
        fh.write("Experiment: %s\n"%self.fnExperiment.get())
        fh.write("Predictor (X): ")
        self.predictor._printToStream(fh)
        fh.write("Predicted (Y): ")
        self.predicted._printToStream(fh)
        fh.write("Analysis: %s\n"%self.analysisDescription)
        fh.write("\n")

        fh.write("[POPULATION PARAMETERS] =============\n")
        fh.write(' '.join(self.analysisParameters)+"\n")
        i=0
        for sampleAnalysis in self.sampleAnalyses:
            outputStr = ""
            j=0
            for parameter in sampleAnalysis.parameters:
                outputStr += "%f "%parameter
                if i==0 and j==0:
                    observations = np.zeros([len(self.sampleAnalyses),len(sampleAnalysis.parameters)])
                observations[i,j]=parameter
                j+=1
            fh.write(outputStr+"\n")
            i+=1
        fh.write("\n")

        mu=np.mean(observations,axis=0)
        C=np.cov(np.transpose(observations))
        sigma = np.sqrt(np.diag(C))
        R=np.corrcoef(np.transpose(observations))
        fh.write("Mean   parameters  = %s\n"%np.array_str(mu))
        fh.write("Median parameters  = %s\n"%np.array_str(np.median(observations,axis=0)))
        fh.write("Lower bound (95%%, independent Gaussians) = %s\n"%np.array_str(mu-1.96*sigma))
        fh.write("Upper bound (95%%, independent Gaussians) = %s\n"%np.array_str(mu+1.96*sigma))
        fh.write("Covariance matrix  =\n%s\n"%np.array_str(C,max_line_width=120))
        fh.write("Correlation matrix  =\n%s\n"%np.array_str(R,max_line_width=120))
        fh.write("\n")

        fh.write("[SAMPLE ANALYSES] ===================\n")
        for sampleAnalysis in self.sampleAnalyses:
            sampleAnalysis._printToStream(fh)

    def load(self,fnAnalysis):
        fh=open(fnAnalysis)
        if not fh:
            raise Exception("Cannot open %s"%fnAnalysis)
        if not verifyMD5(fnAnalysis):
            raise Exception("The file %s has been modified since its creation"%fnAnalysis)
        self.fnAnalysis.set(fnAnalysis)

        for line in fh.readlines():
            line=line.strip()
            if line=="":
                if state==PKPDSignalAnalysis.READING_SAMPLEANALYSIS_PARAMETERS:
                     state=PKPDSignalAnalysis.READING_SAMPLEANALYSIS_NAME
                continue
            if line.startswith('[') and line.endswith('='):
                section = line.split('=')[0].strip().lower()
                if section=="[signal analysis]":
                    state=PKPDSignalAnalysis.READING_EXPERIMENT
                    self.summaryLines.append(line)
                elif section=="[population parameters]":
                    state=PKPDSignalAnalysis.READING_POPULATION_HEADER
                    self.summaryLines.append(line)
                elif section=="[sample analyses]":
                    state=PKPDSignalAnalysis.READING_SAMPLEANALYSIS_NAME
                else:
                    print("Skipping: ",line)

            elif state==PKPDSignalAnalysis.READING_EXPERIMENT:
                tokens = line.split(':')
                self.fnExperiment.set(tokens[1].strip())
                state = PKPDSignalAnalysis.READING_PREDICTOR
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_PREDICTOR:
                tokens = line.split(':')
                self.predictor = PKPDVariable()
                self.predictor.parseTokens(tokens[1].split(';'))
                state = PKPDSignalAnalysis.READING_PREDICTED
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_PREDICTED:
                tokens = line.split(':')
                self.predicted = PKPDVariable()
                self.predicted.parseTokens(tokens[1].split(';'))
                state = PKPDSignalAnalysis.READING_MODEL
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_MODEL:
                tokens = line.split(':')
                self.analysisDescription = tokens[1].strip()
                state = PKPDSignalAnalysis.READING_POPULATION
                self.summaryLines.append(line)
                self.summaryLines.append("\n")

            elif state==PKPDSignalAnalysis.READING_POPULATION_HEADER:
                self.analysisParameters=line.split(' ')
                state = PKPDSignalAnalysis.READING_POPULATION

            elif state==PKPDSignalAnalysis.READING_POPULATION:
                self.summaryLines.append(line)

            elif state==PKPDSignalAnalysis.READING_SAMPLEANALYSIS_NAME:
                tokens = line.split(':')
                self.sampleAnalyses.append(PKPDSampleSignalAnalysis())
                self.sampleAnalyses[-1].sampleName = tokens[1].strip()
                self.sampleAnalyses[-1].parameters=[]
                state = PKPDSignalAnalysis.READING_SAMPLEANALYSIS_PARAMETERS

            elif state==PKPDSignalAnalysis.READING_SAMPLEANALYSIS_PARAMETERS:
                tokens=line.split('=')
                self.sampleAnalyses[-1].analysisVariables.append(tokens[0].strip())
                self.sampleAnalyses[-1].parameters.append(float(tokens[1].strip()))

        fh.close()

    def getSampleAnalysis(self,sampleName):
        for sampleAnalysis in self.sampleAnalyses:
            if sampleAnalysis.sampleName == sampleName:
                return sampleAnalysis
        return None
