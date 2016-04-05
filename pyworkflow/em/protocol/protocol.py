# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
In this module are protocol base classes related to EM.
Them should be sub-classes in the different sub-packages from
each EM-software packages.
"""

from itertools import izip

from pyworkflow.protocol import Protocol
from pyworkflow.object import Set
from pyworkflow.em.constants import RELATION_SOURCE, RELATION_TRANSFORM
from pyworkflow.utils.path import cleanPath
from pyworkflow.mapper.sqlite_db import SqliteDb



class EMProtocol(Protocol):
    """ Base class to all EM protocols.
    It will contains some common functionalities. 
    """
    _base = True
    
    def __init__(self, **kwargs):
        Protocol.__init__(self, **kwargs)
        
    def _defineSourceRelation(self, srcObj, dstObj):
        """ Add a DATASOURCE relation between srcObj and dstObj """
        self._defineRelation(RELATION_SOURCE, srcObj, dstObj)
    
    def _defineTransformRelation(self, srcObj, dstObj):
        self._defineRelation(RELATION_TRANSFORM, srcObj, dstObj)
        # A transform relation allways implies a source relation
        self._defineSourceRelation(srcObj, dstObj)
        
    def _insertChild(self, key, child):
        if isinstance(child, Set):
            child.write()
        Protocol._insertChild(self, key, child)
        
    def __str__(self):
        return self.getObjLabel()
    
    def allowsDelete(self, obj):
        return False
        
    
