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
from pk_models import PK_Twocompartments
import pyworkflow.em.biopharmaceutics

# TESTED in test_workflow_gabrielsson_pk07.py
# TESTED in test_workflow_gabrielsson_pk08.py
# TESTED in test_workflow_gabrielsson_pk09.py
# TESTED in test_workflow_gabrielsson_pk10.py
# TESTED in test_workflow_gabrielsson_pk11.py
# TESTED in test_workflow_gabrielsson_pk12.py
# TESTED in test_workflow_gabrielsson_pk13.py
# TESTED in test_workflow_gabrielsson_pk14.py
# TESTED in test_workflow_gabrielsson_pk16.py
# TESTED in test_workflow_gabrielsson_pk19.py
# TESTED in test_workflow_gabrielsson_pk39.py

class ProtPKPDTwoCompartments(ProtPKPDODEBase):
    """ Fit a two-compartmentx model to a set of measurements (any arbitrary dosing regimen is allowed)\n
        The central compartment is referred to as C, while the peripheral compartment as Cp.
        The differential equation is V dC/dt = -(Cl+Clp) * C + Clp * Cp + dD/dt, Vp dCp/dt = Clp * C - Clp * Cp\n
        where C is the concentration of the central compartment, Cl the clearance, V and Vp the distribution volume of the central and peripheral compartment,
        Clp is the distribution rate between the central and the peripheral compartments, and D the input dosing regime.
        Confidence intervals calculated by this fitting may be pessimistic because it assumes that all model parameters
        are independent, which are not. Use Bootstrap estimates instead.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'pk two-compartments'

    def __init__(self,**kwargs):
        ProtPKPDODEBase.__init__(self,**kwargs)

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        self._defineParams1(form, True, "t", "Cp")
        form.addParam('bounds', params.StringParam, label="Parameter bounds ([tlag], sourceParameters, Cl, V, Clp, Vp)", default="",
                      help="Bounds for time delay, central clearance and volume and peripheral clearance and volume. "\
                      'Make sure that the bounds are expressed in the expected units (estimated from the sample itself).'\
                      'Be careful that Cl bounds must be given here. If you have an estimate of the elimination rate, this is Ke=Cl/V. Consequently, Cl=Ke*V ')

    def createModel(self):
        return PK_Twocompartments()
