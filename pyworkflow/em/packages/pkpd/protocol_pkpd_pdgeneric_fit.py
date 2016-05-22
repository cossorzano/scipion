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

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from protocol_pkpd_fit_base import ProtPKPDFitBase
from pd_models import PDLinear


class ProtPKPDGenericFit(ProtPKPDFitBase):
    """ Fit a generic model. The observed measurement is modelled as Y=f(X).\n
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'fit pd generic'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        self._defineParams1(form,"Cp","E")
        form.addParam('modelType', params.EnumParam, choices=["Linear","Log-linear","Saturated","Sigmoid","Gompertz",
                                                              "Logistic1 ","Logistic 2","Logistic 3","Logistic 4",
                                                              "Richards","Morgan-Mercer-Flodin","Weibull"],
                      label="Generic model", default=0,
                      help='Linear: Y=e0+s*X\nLog-linear: Y=m*log(X-X0)\n')
        form.addParam('fitType', params.EnumParam, choices=["Linear","Logarithmic","Relative"], label="Fit mode", default=1,
                      help='Linear: sum (Cobserved-Cpredicted)^2\nLogarithmic: sum(log10(Cobserved)-log10(Cpredicted))^2\n'\
                           "Relative: sum ((Cobserved-Cpredicted)/Cobserved)^2")
        form.addParam('bounds', params.StringParam, label="Amplitude and time constant bounds", default="", expertLevel=LEVEL_ADVANCED,
                      help='Parameter values for the simulation.\nExample: (1,10);(0,0.05) is (1,10) for the first parameter, (0,0.05) for the second parameter\n'
                           'Linear: e0;s\n'\
                           'Log-linear: m;X0\n')
        form.addParam('confidenceInterval', params.FloatParam, label="Confidence interval=", default=95, expertLevel=LEVEL_ADVANCED,
                      help='Confidence interval for the fitted parameters')
        form.addParam('reportX', params.StringParam, label="Evaluate at X=", default="", expertLevel=LEVEL_ADVANCED,
                      help='Evaluate the model at these X values\nExample 1: [0,5,10,20,40,100]\nExample 2: 0:2:10, from 0 to 10 in steps of 2')

    #--------------------------- STEPS functions --------------------------------------------
    def createModel(self):
        if self.modelType.get()==0:
            return PDLinear()
