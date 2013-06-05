# Enable the 'Prefer executable CLIs' option in Slicer under Modules before running the wizard
from __main__ import qt, ctk, slicer

from CMRToolkitWizardStep import *
from Helper import *
import PythonQt

class CMRToolkitWizardAntrumCutStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '6. PV Antrum Cut' )
    self.setDescription( 'Remove pulmonary veins from endocardial segmentation label image' )

    self.__parent = super( CMRToolkitWizardAntrumCutStep, self )

  def createUserInterface( self ):
    '''
    '''
    self.__layout = self.__parent.createUserInterface()
    print "createUserInterface: CMRToolkitWizardAntrumCutStep"
    
    ##TODO: Allow user to select the wall and endo or should the images be preloaded?
    
    #epiVolumeLabel = qt.QLabel( 'Epi Segmentation (Dilated) Label Image:' )
    #self.__epiVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__epiVolumeSelector.toolTip = "Select the epicardium segmentation image"
    #self.__epiVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__epiVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__epiVolumeSelector.addEnabled = False
    #self.__epiVolumeSelector.renameEnabled = True
    #self.__layout.addRow( epiVolumeLabel, self.__epiVolumeSelector )
    
    #endoVolumeLabel = qt.QLabel( 'Endo Segmentation Label Image:' )
    #self.__endoVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__endoVolumeSelector.toolTip = "Select the endocardium segmentation image"
    #self.__endoVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__endoVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__endoVolumeSelector.addEnabled = False
    #self.__endoVolumeSelector.renameEnabled = True
    #self.__layout.addRow( endoVolumeLabel, self.__endoVolumeSelector )

    antrumCutButton = qt.QPushButton('Load PV Antrum Cut Module')
    self.__layout.addRow(antrumCutButton)
    antrumCutButton.connect('clicked()', self.startAntrumCut)
    
  def startAntrumCut(self):
    # Load epi and endo segmentation volumes from previous steps
    endoSegVolumeID = pNode.GetParameter('endoSegVolumeID')
    endoSegVolume = Helper.getNodeByID(endoSegVolumeID)
    
    wallCleanupSegVolumeID = pNode.GetParameter('wallCleanupSegID')
    wallCleanupSegVolume = Helper.getNodeByID(wallCleanupSegVolumeID)
      
    volumesLogic = slicer.modules.volumes.logic()
    global antrumCutOutputVolume 
    antrumCutOutputVolume = volumesLogic.CreateAndAddLabelVolume( slicer.mrmlScene, endoSegVolume, 'Endocardium No Veins Segmentation Output' )
    
    if endoSegVolume != None and wallCleanupSegVolume != None and endoSegVolumeID != wallCleanupSegVolumeID:
      param = {}
      param['endoLayer'] = endoSegVolume.GetID()
      param['wallLayer'] = wallCleanupSegVolume.GetID()
      param['endoNoVeins'] = antrumCutOutputVolume.GetID()
      slicer.cli.run( slicer.modules.carmapvantrumcut, None, param, wait_for_completion=True )
    else:
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'PV Antrum Cut', 'Both endocardial and wall segmentation images are required to run PV Antrum Cut.')
    
  def validate( self, desiredBranchId ):    
    self.__parent.validate( desiredBranchId )
    
    if antrumCutOutputVolume != None:
      antrumCutVolumeID = antrumCutOutputVolume.GetID()    
      pNode = self.parameterNode()
      pNode.SetParameter('antrumCutVolumeID', antrumCutVolumeID)
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Endocardium No Veins segmentation image not created successfully.')    

  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardAntrumCutStep, self).onEntry(comingFrom, transitionType)
    # Make global in order to be accessible by 'startAntrumCut'
    global pNode 
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

  def onExit(self, goingTo, transitionType):
    pNode = self.parameterNode()
    
    if goingTo.id() != 'EndoNVIsosurface':
      return
      
    super(CMRToolkitWizardAntrumCutStep, self).onExit(goingTo, transitionType)

