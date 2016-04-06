# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules contains basic hierarchy
for EM data objects like: Image, SetOfImage and others
"""

import json
import sys

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

class PKPDExperiment(EMObject):
    READING_GENERAL = 1
    READING_VARIABLES = 2
    READING_DOSES = 3
    READING_SAMPLES = 4
    READING_MEASUREMENTS = 5

    def __init__(self, **args):
        EMObject.__init__(self, **args)
        self.general = {}
        self.variables = {}
        self.samples = {}
        self.dose = {}
        self.measurements = {}

    def load(self, fnExperiment):
        fh=open(fnExperiment,'r')
        if not fh:
            raise Exception("Cannot open the file "+fnExperiment)

        for line in fh.readlines():
            line=line.strip()
            if line=="":
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
                    print("Skipping: ",line)
                    continue
                varname = tokens[0].strip()
                self.variables[varname] = PKPDVariable(tokens)

        fh.close()
        self._printToStream(sys.stdout)

    def write(self, fnExperiment):
        fh=open(fnExperiment,'w')
        self._printToStream(fh)
        fh.close()

    def _printToStream(self,fh):
        fh.write("[EXPERIMENT] ===========================\n")
        for key, value in self.general.iteritems():
            fh.write("%s = %s\n"%(key,value))
        fh.write("\n")

        fh.write("[VARIABLES] ============================\n")
        for key, value in self.variables.iteritems():
            value._printToStream(fh)
        fh.write("\n")
