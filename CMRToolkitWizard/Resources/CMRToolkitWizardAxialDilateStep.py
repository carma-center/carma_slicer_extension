# Enable the 'Prefer executable CLIs' option in Slicer under Modules before running the wizard
from __main__ import qt, ctk, slicer

from CMRToolkitWizardStep import *
from Helper import *
import PythonQt

class CMRToolkitWizardAxialDilateStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '3. Axial dilation of endocardium segmentation' )
    self.setDescription( 'Generate a new epicardium segmentation label image from the endocardium segmentation previously created.' )

    self.__parent = super( CMRToolkitWizardAxialDilateStep, self )

  def createUserInterface( self ):
    '''
    '''
    self.__layout = self.__parent.createUserInterface()
        
    #epiVolumeLabel = qt.QLabel( 'Epi Segmentation (Dilated) Image:' )
    #self.__epiVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__epiVolumeSelector.toolTip = "Create a new output epicardium (dilated) segmentation image"
    #self.__epiVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__epiVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__epiVolumeSelector.addEnabled = True
    #self.__epiVolumeSelector.renameEnabled = True
    #self.__epiVolumeSelector.baseName = "Epi Segmentation"
    #self.__layout.addRow( epiVolumeLabel, self.__epiVolumeSelector )

    #self.__epiVolumeSelector.setCurrentNode( epiOutputVolume )
    
    #endoSegVolumeID = pNode.GetParameter('endoSegVolumeID')
    #endoSegVolume = Helper.getNodeByID(endoSegVolumeID)
    
    axialDilateButton = qt.QPushButton('Load Axial Dilate Module')
    self.__layout.addRow(axialDilateButton)
    axialDilateButton.connect('clicked()', self.startAxialDilate)

  def startAxialDilate(self):
    # Load endo segmentation volume from previous step
    endoSegVolumeID = pNode.GetParameter('endoSegVolumeID')
    endoSegVolume = Helper.getNodeByID(endoSegVolumeID)
    
    volumesLogic = slicer.modules.volumes.logic()
    global epiOutputVolume 
    epiOutputVolume = volumesLogic.CreateAndAddLabelVolume( slicer.mrmlScene, endoSegVolume, 'Epi Segmentation Output' )

    # Initialize the input volumes to run the Axial Dilate Module
    #epiOutputVolume = self.__epiVolumeSelector.currentNode()
      
    if epiOutputVolume != None and endoSegVolume != None:
      param = {}
      param['targetFileName'] = endoSegVolume.GetID()
      param['outputFileName'] = epiOutputVolume.GetID()
      slicer.cli.run( slicer.modules. cmrtoolkitaxialdilate, None, param, wait_for_completion=True )
    else:
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Axial Dilate', 'Output epicardial segmentation image is required to run Axial Dilate.')
    
  def validate( self, desiredBranchId ):    
    self.__parent.validate( desiredBranchId )
    #epiOutputVolume = self.__epiVolumeSelector.currentNode()
    
    if epiOutputVolume != None:
      epiSegID = epiOutputVolume.GetID()    
      pNode = self.parameterNode()
      pNode.SetParameter('epiSegVolumeID', epiSegID)        
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Please select an epi segmentation image to proceed.')    

  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardAxialDilateStep, self).onEntry(comingFrom, transitionType)
    # Make pNode global in order to be accessible by 'startAxialDilate'
    global pNode 
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)
    
  def onExit(self, goingTo, transitionType):
    pNode = self.parameterNode()
    
    #if goingTo.id() != 'BooleanRemove':
    #  return
      
    super(CMRToolkitWizardAxialDilateStep, self).onExit(goingTo, transitionType)

      