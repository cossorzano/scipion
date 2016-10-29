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

import math
import sys
import numpy as np

import pyworkflow.protocol.params as params
from pyworkflow.em.protocol.protocol_pkpd import ProtPKPD
from pyworkflow.em.data import PKPDVariable

class ProtPKPDSimulateDrugInteractions(ProtPKPD):
    """ Simulate drug interactions as recommended in EMA CHMP/EWP/560/95 \n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'simulate drug interactions'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form, fullForm=True):
        form.addSection('Input')
        fromTo = form.addLine('Inhibitor range [I]',
                           help='[I] is the maximal total (free and bound) systemic inhibitor concentration in plasma at the highest dose')
        fromTo.addParam('I0', params.FloatParam, default=0, label='Min (uM)')
        fromTo.addParam('IF', params.FloatParam, default=10, label='Max (uM)')

        form.addParam('Ki', params.StringParam, default="5", label='Ki (uM)',
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")
        form.addParam("doBasic", params.BooleanParam, default=True, label="Basic model",
                      help="Investigational drug likely to be a reversible inhibitor if [I]/Ki>0.02")


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate')

    #--------------------------- STEPS functions --------------------------------------------
    def parseList(self, strList):
        return [float(v) for v in strList.split(' ')]

    def runSimulate(self):
        KiList = self.parseList(self.Ki.get())
        I = np.arange(self.I0.get(), self.IF.get(), (self.IF.get()-self.I0.get())/100)
        if self.doBasic:
            R1basic = np.empty((len(KiList),I.size))
            for i,Ki in enumerate(KiList):
                R1basic[i,:] = 1+I/Ki

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        return msg

    def _validate(self):
        msg=[]
        return msg

    def _citations(self):
        return ['CHMPEWP56095']