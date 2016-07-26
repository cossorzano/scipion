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
# *  All comments concerning this program package may be sent to thes
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************

import pyworkflow.protocol.params as params
from protocol_pkpd_ode_base import ProtPKPDODEBase
from pk_models import PK_Monocompartment
import biopharmaceutics


class ProtPKPDEV1MonoCompartment(ProtPKPDODEBase):
    """ Fit a monocompartmental model with first order absorption to a set of measurements obtained by oral doses (any arbitrary dosing regimen is allowed)\n
        The differential equation is D(t)=Amax*(1-exp(-Ka*t)) and dC/dt = -Cl * C + 1/V * dD/dt\n
        where C is the concentration, Cl the clearance, V the distribution volume, and D the input dosing regime. Ka is the absorbption rate and
        Amax the dose (maximum amount).
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'ev1 monocompartment'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form, True, "t", "Cp")
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], Ka, Cl, V)", default="",
                      help="Bounds for the tlag (if it must be estimated), absorption, clearance and volume. Example: (0.01,0.04);(0.2,0.4);(10,20). "\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).'\
                      'If tlag must be estimated, its bounds must always be specified. '\
                      'Be careful that Cl bounds must be given here. If you have an estimate of the elimination rate, this is Ke=Cl/V. Consequently, Cl=Ke*V ')

    def configureSource(self):
        self.drugSource.type = biopharmaceutics.DrugSource.EV
        self.drugSource.evProfile = biopharmaceutics.BiopharmaceuticsModelOrder1()
        self.drugSource.evProfile.Amax = self.drugSource.getCumulatedDose(self.model.t0,self.model.tF)

    def createModel(self):
        return PK_Monocompartment()

