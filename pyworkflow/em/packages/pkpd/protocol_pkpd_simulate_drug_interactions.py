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
                           help='[I] is [I]gut = Molar Dose/250mL.')
        fromTo.addParam('I0', params.FloatParam, default=0, label='Min (uM)')
        fromTo.addParam('IF', params.FloatParam, default=10, label='Max (uM)')

        form.addParam("doReversible", params.BooleanParam, default=False, label="Reversible inhibition",
                      help="Investigational drug likely to be a reversible inhibitor if R1=1+[I]/Ki>1.02")
        form.addParam('KiReversible', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doReversible",
                      help="Ki is the in vitro unbound reversible inhibition constant. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doTimeDependent", params.BooleanParam, default=False, label="Time dependent inhibition",
                      help="Investigational drug likely to be a time dependent inhibitor if R2=1+kinact/kdeg*[I]/([I]+Ki)>1.25")
        form.addParam("kinact",params.StringParam, default="0.01", label='Max. inactivation rate (kinact [min^-1])', condition="doTimeDependent",
                      help="Maximal inactivation rate. Several constants can be given separated by space, e.g., 0.01 0.005")
        form.addParam("kdeg",params.StringParam, default="0.008", label='Apparent first order degradation rate (kdeg [min^-1])', condition="doTimeDependent",
                      help="kdeg is the apparent first order degradation rate constant of the affected enzyme. Several constants can be given separated by space, e.g., 0.01 0.005")
        form.addParam('KiTime', params.StringParam, default="5", label='Inhibition constant (Ki [uM])', condition="doTimeDependent",
                      help="KI is the inhibitor concentration which yields 50% of the maximum inactivation rate. Several constants can be given separated by space, e.g., 5 10")

        form.addParam("doInduction", params.BooleanParam, default=False, label="Induction",
                      help="Investigational drug likely to be a time dependent inhibitor if R2=1+kinact/kdeg*[I]/([I]+Ki)>1.25")
        form.addParam("d",params.StringParam, default="1", label='Scaling factor', condition="doInduction",
                      help="Several constants can be given separated by space, e.g., 1 1.05")
        form.addParam("Emax",params.StringParam, default="1", label='Max. Induction effect (Emax)', condition="doInduction",
                      help="Several constants can be given separated by space, e.g., 1 2")
        form.addParam('EC50', params.StringParam, default="5", label='Half Max. Effect Conc (EC50 [uM])', condition="doInduction",
                      help="EC50 is the concentration causing half maximal effect. Several constants can be given separated by space, e.g., 5 6")

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSimulate')

    #--------------------------- STEPS functions --------------------------------------------
    def parseList(self, strList):
        return [float(v) for v in strList.split(' ')]

    def runSimulate(self):
        I = np.arange(self.I0.get(), self.IF.get(), (self.IF.get()-self.I0.get())/100)
        R = []
        Rlegends = []
        if self.doReversible:
            KiList = self.parseList(self.KiReversible.get())
            R1 = np.empty((len(KiList),I.size))
            for Ki in KiList:
                legend="Rev. Inh. Ki=%f [uM]"%Ki
                print("Simulating %s"%legend)
                R.append(1+I/Ki)
                Rlegends.append(legend)

        if self.doTimeDependent:
            KiList = self.parseList(self.KiReversible.get())
            kdegList = self.parseList(self.kdeg.get())
            kinactList = self.parseList(self.kinact.get())
            for Ki in KiList:
                for kdeg in kdegList:
                    for kinact in kinactList:
                        legend="Time Dep. Inh. Ki=%f [uM], kdeg=%f [min^-1], kinact=%f [min^-1]"%(Ki,kdeg,kinact)
                        print("Simulating %s"%legend)
                        R.append(1+kinact/kdeg*I/(Ki+I))
                        Rlegends.append(legend)

        if self.doInduction:
            EC50List = self.parseList(self.EC50.get())
            EmaxList = self.parseList(self.Emax.get())
            dList = self.parseList(self.d.get())
            for EC50 in EC50List:
                for Emax in EmaxList:
                    for d in dList:
                        legend="Induction EC50=%f [uM], Emax=%f, d=%f"%(EC50,Emax,d)
                        print("Simulating %s"%legend)
                        R.append(1/(1+d*Emax*I/(EC50+I)))
                        Rlegends.append(legend)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        msg=[]
        return msg

    def _validate(self):
        msg=[]
        return msg

    def _citations(self):
        return ['CHMPEWP56095']