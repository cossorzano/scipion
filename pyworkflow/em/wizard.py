# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
This module implement some wizards
"""

import os
from os.path import basename, exists
import Tkinter as tk
import ttk

from pyworkflow.wizard import Wizard
import pyworkflow.gui.dialog as dialog
from pyworkflow.gui.widgets import LabelSlider
from pyworkflow.gui.tree import BoundTree, TreeProvider
from pyworkflow import findResource
from pyworkflow.object import PointerList

from pyworkflow.em.convert import ImageHandler


#===============================================================================
#    Wizard EM base class
#===============================================================================

class EmWizard(Wizard):    
    
    def _getText(self, obj):
        index = obj.getIndex()
        text = os.path.basename(obj.getFileName())
        if index:
            return "%03d@%s" % (index, text)
        return text
    
    
    def _getListProvider(self, objs):
        """ This should be implemented to return the list
        of object to be displayed in the tree.
        """
        provider = None
        if objs.hasValue():
            # If objs is a PointerList currently it can only be formed of SetOfVolumes and Volume
            # (for protocol align_volume). Should this change review this part
            if isinstance(objs, PointerList):
                vols_total = []
                for pointer in objs:
                    obj = pointer.get()
                    print obj
                    vols = self._getVols(obj)
                    vols_total.extend(vols)
                provider = ListTreeProvider(vols_total)
            else:
                objs = objs.get()

                if isinstance(objs, SetOfMicrographs):
                    mics = self._getMics(objs)
                    provider = ListTreeProvider(mics)

                if isinstance(objs, SetOfParticles):
                    particles = self._getParticles(objs)
                    provider = ListTreeProvider(particles)
                    provider.getText = self._getText

                if isinstance(objs, SetOfVolumes) or isinstance(objs, Volume):
                    vols = self._getVols(objs)
                    provider = ListTreeProvider(vols)
            
            return provider
        return None
    
    def _getInputProtocol(self, targets, protocol):
        label = []
        value = []
        
        for k, v in targets:
            if k.__name__ == protocol.getClassName():
                label = v
                for val in v:
                    value.append(protocol.getAttributeValue(val))
        
        if len(label) > 1:     
            return label, value
        else:
            return label[0], value[0]
    
#===============================================================================
#    Provider class (for objects used in the wizards)
#===============================================================================

class ListTreeProvider(TreeProvider):
    """ Simple list tree provider. """
    def __init__(self, objList=None):
        self.objList = objList
        self.getColumns = lambda: [('Object', 150)]
        self.getObjects = lambda: self.objList
    
    def getObjectInfo(self, obj):
        info = {'key': obj.getObjId(), 'text': self.getText(obj), 'values': ()}
        return info
    
    def getText(self, obj):
        """ Get the text to display for an object. """
        return os.path.basename(obj.getFileName())
    
    def getObjs(self):
        """ Get the objects. """
        return self.objList


