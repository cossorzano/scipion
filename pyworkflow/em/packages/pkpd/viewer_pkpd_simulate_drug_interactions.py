# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com)
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

from itertools import izip
import numpy as np

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.plotter import EmPlotter

from protocol_pkpd_simulate_drug_interactions import ProtPKPDSimulateDrugInteractions

class PKPDSimulateDrugInteractionsViewer(Viewer):
    _targets = [ProtPKPDSimulateDrugInteractions]
    _environments = [DESKTOP_TKINTER]

    def visualize(self, obj, **kwargs):
        prot = obj
        fnProfiles = prot._getPath("profiles.txt")
        fh = open(fnProfiles,"r")
        Rlegends = []
        R = []
        state = 0
        for line in fh:
            if state==0:
                Rlegends.append(line.strip())
                Ri=[]
                state=1
            elif state==1:
                tokens=line.strip().split()
                if len(tokens)==0:
                    R.append(Ri)
                    state=0
                else:
                    if len(Ri)==0:
                        for n in range(len(tokens)):
                            Ri.append([])
                for n in range(len(tokens)):
                    Ri[n].append(tokens[n])
        fh.close()

        plotter = None
        for legend, Ri  in izip(Rlegends, R):
            if plotter is None:
                plotter = EmPlotter()
                doShow = True
                ax = plotter.createSubPlot("Plot", "[I]", "R")
            else:
                doShow = False
                ax = plotter.getLastSubPlot()

            x = np.asarray(Ri[0])
            if len(Ri)==2:
                y=np.asarray(Ri[1])
            else:
                y=np.asarray(Ri[2])
            ax.plot(x, y, label=legend)
            leg = ax.legend()
            if leg:
                leg.draggable()

            if doShow:
                plotter.show()
            else:
                plotter.draw()
