# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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
This module implement viewers for some type of common objects.
"""
from __future__ import print_function
import os
import sys
import shlex
import ast
from threading import Thread
from multiprocessing.connection import Client
from numpy import flipud
import socket

from pyworkflow.viewer import View, Viewer, CommandView, DESKTOP_TKINTER
from pyworkflow.utils import Environ, runJob
from pyworkflow.utils import getFreePort
from pyworkflow.gui.matplotlib_image import ImageWindow

# From pyworkflow.em level
import showj
import xmipp

#------------------------ Some common Views ------------------

class DataView(View):
    """ Wrapper the arguments to showj (either web or desktop). """
    def __init__(self, path, viewParams={}, **kwargs):
        View.__init__(self)
        self._memory = '2g'
        self._loadPath(path)
        self._env = kwargs.get('env', {})
        self._viewParams = viewParams
            
    def _loadPath(self, path):
        # Check if there is a table name with @ in path
        # in that case split table name and path
        # table names can never starts with a number
        # this is considering an image inside an stack
        if '@' in path and path[0] not in '0123456789':
            self._tableName, self._path = path.split('@')
        else:
            self._tableName, self._path = None, path
            
    def show(self):        
        showj.runJavaIJapp(self._memory, 'xmipp.viewer.scipion.ScipionViewer',
                           self.getShowJParams(), env=self._env)
    
    def getShowJParams(self):
        tableName = '%s@' % self._tableName if self._tableName else ''
        params = '-i "%s%s"' % (tableName, self._path)
        for key, value in self._viewParams.items():
            params = "%s --%s %s"%(params, key, value)
        
        return params
    
    def getShowJWebParams(self):
    
    #=OLD SHOWJ WEB DOCUMENTATION===============================================
    # Extra parameters can be used to configure table layout and set render function for a column
    # Default layout configuration is set in ColumnLayoutProperties method in layout_configuration.py
    # 
    # Parameters are formed by: [label]___[property]: [value]. E.g.: id___visible:True or micrograph___renderFunc:"get_image_psd"
    # Properties to be configured are:
    #    visible: Defines if this column is displayed
    #    allowSetVisible: Defines if user can change visible property (show/hide this column).
    #    editable: Defines if this column is editable, ie user can change field value.
    #    allowSetEditable: Defines if user can change editable property (allow editing this column).
    #    renderable: Defines if this column is renderizable, ie it renders data column using renderFunc
    #    allowSetRenderable: Defines if user can change renderable property.
    #    renderFunc: Function to be used when this field is rendered. (it has to be inserted in render_column method)
    #    extraRenderFunc: Any extra parameters needed for rendering. Parameters are passed like in a url ie downsample=2&lowPass=3.5
    # 
    # Example:
    # extraParameters["id___visible"]=True
    # extraParameters["micrograph___renderFunc"]="get_image_psd"
    # extraParameters["micrograph___extraRenderFunc"]="downsample=2"
    #===========================================================================
    
        parameters = {
            showj.MODE, # FOR MODE TABLE OR GALLERY
            showj.VISIBLE,
            showj.ZOOM,
            showj.ORDER,
            showj.RENDER,
            showj.SORT_BY
#             'columns',
        }
        
        params = {}
        
        for key, value in self._viewParams.items():
            print (str(key), ":",str(value))
            if key in parameters:
                if key == 'mode' and value =='metadata':
                    value = 'table'
                params[key] = value
        
        return params
        
    def getPath(self):
        return self._path
    
    def getTableName(self):
        return self._tableName
        
        
class ObjectView(DataView):
    """ Wrapper to DataView but for displaying Scipion objects. """
    def __init__(self, project, inputid, path, other='', viewParams={}, **kwargs):
        DataView.__init__(self, path, viewParams, **kwargs)
        self.type = type
        self.port = project.port
        self.inputid = inputid
        self.other = other
        
    def getShowJParams(self):
        # mandatory to provide scipion params
        params = DataView.getShowJParams(self) + ' --scipion %s %s %s' % (self.port, self.inputid, self.other)
        return params
    
    def show(self):
        showj.runJavaIJapp(self._memory, 'xmipp.viewer.scipion.ScipionViewer', self.getShowJParams(), env=self._env)
        
