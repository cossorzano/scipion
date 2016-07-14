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

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDExperiment, PKPDVariable, PKPDDose, PKPDSample
import sys


class ProtPKPDImportFromCSV(ProtPKPD):
    """ Import experiment from CSV.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'import from csv'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputFile', params.PathParam,
                      label="File path", allowsNull=False,
                      help='Specify a path to desired CSV file.')
        form.addParam('title', params.StringParam, label="Title", default="My experiment")
        form.addParam('comment', params.StringParam, label="Comment", default="")
        form.addParam('variables', params.TextParam, label="Variables", default="",
                      help="Structure: [Variable Name] ; [Units] ; [Type] ; [Comment]\n"\
                           "The variable name should have no space or special character\n"\
                           "Valid units are: h, mg, ug, ug/mL, ...\n"\
                           "Type is either numeric or text\n"\
                           "The comment may be empty\n"\
                           "\nIt is important that there are three semicolons (none of them may be missing even if the comment is not present).\n"\
                           "Examples:\n"\
                           "t ; h ; numeric ; time ; \n"\
                           "Cp ; ug/mL ; numeric ; measurement ; plasma concentration\n"\
                           "weight ; g ; numeric; label ; weight of the animal\n"\
                           "sex ; g ; text ; label ; sex of the animal\n")
        form.addParam('doses', params.TextParam, label="Doses", default="",
                      help="Structure: [Dose Name] ; [Description] ; [Units] \n"\
                           "The dose name should have no space or special character\n"\
                           "Valid units are: h, mg, ug, ug/mL, ...\n"\
                           "The description is either a bolus or an infusion as shown in the examples\n"\
                           "\nIt is important that there are two semicolons.\n"\
                           "Examples:\n"\
                           "Infusion0 ; infusion t=0.500000...0.750000 d=60*weight/1000; mg\n"\
                           "Bolus1 ; bolus t=2.000000 d=100; mg\n"\
                           "Bolus0 ; bolus t=0.000000 d=60*weight/1000; mg\n")
        form.addParam('dosesToSamples', params.TextParam, label="Assign doses to samples", default="",
                      help="Structure: [Sample Name] ; [DoseName1,DoseName2,...] \n"\
                           "The sample name should have no space or special character\n"\
                           "\nIt is important that there is one semicolon.\n"\
                           "Examples:\n"\
                           "FemaleRat1 ; Bolus0,Bolus1,Infusion0\n")


    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',self.inputFile.get())

    #--------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self, objId):
        experiment = PKPDExperiment()
        experiment.general["title"]=self.title.get()
        experiment.general["comment"]=self.comment.get()

        # Read the variables
        for line in self.variables.get().split(';;'):
            tokens = line.split(';')
            if len(tokens)!=5:
                print("Skipping variable: ",line)
                continue
            varname = tokens[0].strip()
            experiment.variables[varname] = PKPDVariable()
            experiment.variables[varname].parseTokens(tokens)

        # Read the doses
        for line in self.doses.get().split(';;'):
            tokens = line.split(';')
            if len(tokens)!=3:
                print("Skipping dose: ",line)
                continue
            dosename = tokens[0].strip()
            experiment.doses[dosename] = PKPDDose()
            experiment.doses[dosename].parseTokens(tokens)

        # Read the sample doses
        for line in self.dosesToSamples.get().split(';;'):
            tokens = line.split(';')
            if len(tokens)<2:
                print("Skipping sample: ",line)
                continue
            samplename = tokens[0].strip()
            tokens[1]="dose="+tokens[1]
            experiment.samples[samplename] = PKPDSample()
            experiment.samples[samplename].parseTokens(tokens,experiment.variables, experiment.doses)

        # Read the measurements
        fh=open(self.inputFile.get())
        lineNo = 1
        for line in fh.readlines():
            tokens = line.split(';')
            if len(tokens)==0:
                continue
            if lineNo==1:
                listOfVariables=[]
                listOfSkips=[]
                iSampleName=-1
                varNo=0
                for token in tokens:
                    varName=token.strip()
                    if varName=="SampleName":
                        iSampleName=varNo
                    listOfVariables.append(varName)
                    listOfSkips.append(not (varName in experiment.variables))
                    varNo+=1
                if iSampleName==-1:
                    raise Exception("Cannot find the SampleName in: %s\n"%line)
                print(listOfVariables)
            else:
                if len(tokens)!=len(listOfSkips):
                    print("Skipping line: %s"%line)
                    print("   It does not have the same number of values as the header")
                varNo = 0
                sampleName = tokens[iSampleName].strip()
                if not sampleName in experiment.samples:
                    print("Skipping sample: The sample %s does not have a dose"%varName)
                    continue
                samplePtr=experiment.samples[sampleName]
                for skip in listOfSkips:
                    if not skip:
                        varName = listOfVariables[varNo]
                        varRole = experiment.variables[listOfVariables[varNo]].role
                        if varRole == PKPDVariable.ROLE_LABEL:
                            samplePtr.descriptors[listOfVariables[varNo]] = tokens[varNo]
                        else:
                            ok=True
                            if tokens[varNo]=="NA":
                                ok = (varRole != PKPDVariable.ROLE_TIME)
                            if ok:
                                if not hasattr(samplePtr,"measurement_%s"%varName):
                                    setattr(samplePtr,"measurement_%s"%varName,[])
                                    samplePtr.measurementPattern.append(varName)
                                    print(samplePtr.measurementPattern)
                                exec("samplePtr.measurement_%s.append('%s')"%(varName,tokens[varNo].strip()))
                            else:
                                raise Exception("Time measurements cannot be NA")
                    varNo+=1
            lineNo+=1

        experiment.write(self._getPath("experiment.pkpd"))
        experiment._printToStream(sys.stdout)
        self._defineOutputs(outputExperiment=experiment)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        return ["Input file: %s"%self.inputFile.get()]