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
"""
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""

from os.path import basename, join, exists
import numpy as np

from pyworkflow.utils.path import cleanPath, makePath, cleanPattern
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.utils.process import runJob
from pyworkflow.em.data import PKPDExperiment

from tk_experiment import ExperimentWindow


# from protocol_nma_dimred import XmippProtDimredNMA, DIMRED_MAPPINGS
# from data import Point, Data
# from plotter import XmippNmaPlotter
# from gui import ClusteringWindow, TrajectoriesWindow


class PKPDExperimentViewer(Viewer):
    """ Visualization of a given PKPDExperiment
    """
    _targets = [PKPDExperiment]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def _visualize(self, obj, **kwargs):
        obj.load()
        self.clusterWindow = self.tkWindow(ExperimentWindow,
                                           title='Experiment Viewer',
                                           experiment=obj,
                                           callback=self._createExperiment
                                           )
        return [self.clusterWindow]

    def _createExperiment(self):
        """ Create a new experiment after manipulation of
        the currently displayed experiment and register the action.
        """
        # Write the particles
        prot = self.protocol
        project = prot.getProject()
        inputSet = prot.getInputParticles()
        fnSqlite = prot._getTmpPath('cluster_particles.sqlite')
        cleanPath(fnSqlite)
        partSet = SetOfParticles(filename=fnSqlite)
        partSet.copyInfo(inputSet)
        for point in self.getData():
            if point.getState() == Point.SELECTED:
                particle = inputSet[point.getId()]
                partSet.append(particle)
        partSet.write()
        partSet.close()

        from protocol_batch_cluster import BatchProtNMACluster
        newProt = project.newProtocol(BatchProtNMACluster)
        clusterName = self.clusterWindow.getClusterName()
        if clusterName:
            newProt.setObjLabel(clusterName)
        newProt.inputNmaDimred.set(prot)
        newProt.sqliteFile.set(fnSqlite)

        project.launchProtocol(newProt)

    def _loadAnimationData(self, obj):
        prot = self.protocol
        animationName = obj.getFileName() # assumes that obj.getFileName is the folder of animation
        animationPath = prot._getExtraPath(animationName)
        #animationName = animationPath.split('animation_')[-1]
        animationRoot = join(animationPath, animationName)

        animationSuffixes = ['.vmd', '.pdb', 'trajectory.txt']
        for s in animationSuffixes:
            f = animationRoot + s
            if not exists(f):
                self.errorMessage('Animation file "%s" not found. ' % f)
                return

        # Load animation trajectory points
        trajectoryPoints = np.loadtxt(animationRoot + 'trajectory.txt')
        data = PathData(dim=trajectoryPoints.shape[1])

        for i, row in enumerate(trajectoryPoints):
            data.addPoint(Point(pointId=i+1, data=list(row), weight=1))

        self.trajectoriesWindow.setPathData(data)
        self.trajectoriesWindow.setAnimationName(animationName)
        self.trajectoriesWindow._onUpdateClick()

        def _showVmd():
            vmdFn = animationRoot + '.vmd'
            VmdView(' -e %s' % vmdFn).show()

        self.getTkRoot().after(500, _showVmd)


    def _loadAnimation(self):
        prot = self.protocol
        browser = FileBrowserWindow("Select the animation folder (animation_NAME)",
                                    self.getWindow(), prot._getExtraPath(),
                                    onSelect=self._loadAnimationData)
        browser.show()

    def _generateAnimation(self):
        prot = self.protocol
        projectorFile = prot.getProjectorFile()

        animation = self.trajectoriesWindow.getAnimationName()
        animationPath = prot._getExtraPath('animation_%s' % animation)

        cleanPath(animationPath)
        makePath(animationPath)
        animationRoot = join(animationPath, 'animation_%s' % animation)

        trajectoryPoints = np.array([p.getData() for p in self.trajectoriesWindow.pathData])

        if projectorFile:
            M = np.loadtxt(projectorFile)
            deformations = np.dot(trajectoryPoints, np.linalg.pinv(M))
            np.savetxt(animationRoot + 'trajectory.txt', trajectoryPoints)
        else:
            Y = np.loadtxt(prot.getOutputMatrixFile())
            X = np.loadtxt(prot.getDeformationFile())
            # Find closest points in deformations
            deformations = [X[np.argmin(np.sum((Y - p)**2, axis=1))] for p in trajectoryPoints]

        pdb = prot.getInputPdb()
        pdbFile = pdb.getFileName()

        modesFn = prot.inputNMA.get()._getExtraPath('modes.xmd')

        for i, d in enumerate(deformations):
            atomsFn = animationRoot + 'atomsDeformed_%02d.pdb' % (i+1)
            cmd = '-o %s --pdb %s --nma %s --deformations %s' % (atomsFn, pdbFile, modesFn, str(d)[1:-1])
            runJob(None, 'xmipp_pdb_nma_deform', cmd, env=prot._getEnviron())

        # Join all deformations in a single pdb
        # iterating going up and down through all points
        # 1 2 3 ... n-2 n-1 n n-1 n-2 ... 3, 2
        n = len(deformations)
        r1 = range(1, n+1)
        r2 = range(2, n) # Skip 1 at the end
        r2.reverse()
        loop = r1 + r2

        trajFn = animationRoot + '.pdb'
        trajFile = open(trajFn, 'w')

        for i in loop:
            atomsFn = animationRoot + 'atomsDeformed_%02d.pdb' % i
            atomsFile = open(atomsFn)
            for line in atomsFile:
                trajFile.write(line)
            trajFile.write('TER\nENDMDL\n')
            atomsFile.close()

        trajFile.close()
        # Delete temporary atom files
        cleanPattern(animationRoot + 'atomsDeformed_??.pdb')

        # Generate the vmd script
        vmdFn = animationRoot + '.vmd'
        vmdFile = open(vmdFn, 'w')
        vmdFile.write("""
        mol new %s
        animate style Loop
        display projection Orthographic
        mol modcolor 0 0 Index
        mol modstyle 0 0 Beads 1.000000 8.000000
        animate speed 0.5
        animate forward
        """ % trajFn)
        vmdFile.close()

        VmdView(' -e ' + vmdFn).show()


    def loadData(self):
        """ Iterate over the images and the output matrix txt file
        and create a Data object with theirs Points.
        """
        matrix = np.loadtxt(self.protocol.getOutputMatrixFile())
        particles = self.protocol.getInputParticles()

        data = Data()
        for i, particle in enumerate(particles):
            data.addPoint(Point(pointId=particle.getObjId(),
                                data=matrix[i, :],
                                weight=particle._xmipp_cost.get()))

        return data

