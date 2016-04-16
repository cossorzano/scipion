# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
In this module are protocol base classes related to PKPD

"""
from pyworkflow.em.protocol import *
from pyworkflow.em.data import PKPDExperiment

class ProtPKPD(EMProtocol):
    def printSection(self, msg):
        print("**********************************************************************************************")
        print(msg)
        print("**********************************************************************************************")

    def readExperiment(self):
        experiment = PKPDExperiment()
        experiment.load(self.inputExperiment.get().fnPKPD)
        self.printSection("Reading %s"%self.inputExperiment.get().fnPKPD)
        experiment._printToStream(sys.stdout)
        return experiment

    def writeExperiment(self, experiment, fnOut):
        self.printSection("Writing %s"%fnOut)
        experiment._printToStream(sys.stdout)
        experiment.write(fnOut)
