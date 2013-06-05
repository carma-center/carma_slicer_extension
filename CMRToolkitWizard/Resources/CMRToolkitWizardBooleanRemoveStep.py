# Enable the 'Prefer executable CLIs' option in Slicer under Modules before running the wizard
from __main__ import qt, ctk, slicer

from CMRToolkitWizardStep import *
from Helper import *
import PythonQt

class CMRToolkitWizardBooleanRemoveStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '4. Generate Wall Segmentation' )
    self.setDescription( 'Subtract the endocardial segmentation from the epicardial segmentation to generate the wall segmentation' )

    self.__parent = super( CMRToolkitWizardBooleanRemoveStep, self )

  def createUserInterface( self ):
    '''
    '''
    self.__layout = self.__parent.createUserInterface()
    
    ##TODO: Allow user to select the epi and endo or should the images be preloaded?
    
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

    booleanRemoveButton = qt.QPushButton('Load Boolean Remove Module')
    self.__layout.addRow(booleanRemoveButton)
    booleanRemoveButton.connect('clicked()', self.startBooleanRemove)

  def startBooleanRemove(self):
    # Load epi and endo segmentation volumes from previous steps
    endoSegVolumeID = pNode.GetParameter('endoSegVolumeID')
    endoSegVolume = Helper.getNodeByID(endoSegVolumeID)
    
    epiSegVolumeID = pNode.GetParameter('epiSegVolumeID')
    epiSegVolume = Helper.getNodeByID(epiSegVolumeID)
      
    volumesLogic = slicer.modules.volumes.logic()
    global wallOutputVolume 
    wallOutputVolume = volumesLogic.CreateAndAddLabelVolume( slicer.mrmlScene, epiSegVolume, 'Wall Segmentation Output' )
    
    if epiSegVolume != None and endoSegVolume != None and epiSegVolumeID != endoSegVolumeID:
      param = {}
      param['inputVolume1'] = epiSegVolume.GetID()
      param['inputVolume2'] = endoSegVolume.GetID()
      param['outputVolume'] = wallOutputVolume.GetID()
      slicer.cli.run( slicer.modules.carmabooleanremovefilter, None, param, wait_for_completion=True )
    else:
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Boolean Remove', 'Both epicardial and endocardial segmentation images are required to run Boolean Remove.')
                              
  def validate( self, desiredBranchId ):    
    self.__parent.validate( desiredBranchId )

    if wallOutputVolume != None:
      wallSegID = wallOutputVolume.GetID()    
      pNode = self.parameterNode()
      pNode.SetParameter('wallSegVolumeID', wallSegID)        
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Wall segmentation image not created successfully.')    

  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardBooleanRemoveStep, self).onEntry(comingFrom, transitionType)
    # Make global in order to be accessible by 'startBooleanRemove'
    global pNode 
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)
    
  def onExit(self, goingTo, transitionType):
    pNode = self.parameterNode()

    if goingTo.id() != 'LAWallCleanup':
      return
      
    super(CMRToolkitWizardBooleanRemoveStep, self).onExit(goingTo, transitionType)
    
    
    