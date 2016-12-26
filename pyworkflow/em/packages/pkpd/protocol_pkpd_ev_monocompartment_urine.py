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
from pk_models import PK_MonocompartmentUrine
import pyworkflow.em.biopharmaceutics


class ProtPKPDEV1MonoCompartmentUrine(ProtPKPDODEBase):
    """ Fit a monocompartmental model to a set of measurements obtained by intravenous doses (any arbitrary dosing regimen is allowed)\n
        and urine measurements of cumulated amount of drug
        The differential equation is dC/dt = -Cl * C/V + 1/V * dD/dt and dA/dt = fe * Cl * C\n
        where C is the concentration, Cl the clearance, V the distribution volume, D the input dosing regime and fe the fraction excreted.
Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'ev1 monocompartment urine'

    def __init__(self, **kwargs):
        ProtPKPDODEBase.__init__(self, **kwargs)
        self.ncaProt = None

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form, True, "t", "Cp")
        form.addParam('Au', params.StringParam, label="Cumulated amount in urine", default="Au")
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], Ka, Cl, V, fe)", default="",
                      help="Bounds for the tlag (if it must be estimated), 1st order input rate, clearance, volume and fraction excreted."\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).'\
                      'If tlag must be estimated, its bounds must always be specified')

    def getXYvars(self):
        self.varNameX=self.predictor.get()
        self.varNameY=[self.predicted.get(),self.Au.get()]

    def configureSource(self, drugSource):
        drugSource.type = biopharmaceutics.DrugSource.EV
        drugSource.evProfile = biopharmaceutics.BiopharmaceuticsModelOrder1()

    def createModel(self):
        return PK_MonocompartmentUrine()
