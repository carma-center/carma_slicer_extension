from __main__ import qt, ctk

from LASegmentationWorkflowStep import *
from Helper import *

class LASegmentationWorkflowSelectDataStep( LASegmentationWorkflowStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '1. Select input data' )
    self.setDescription( 'Select Cardiac MRI data to be analyzed. Use the test data or your own images.' )

    self.__parent = super( LASegmentationWorkflowSelectDataStep, self )

  def killButton(self):
    # Hide unneccesary button
    bl = slicer.util.findChildren(text='LAEndo*')
    if len(bl):
      bl[0].hide()

  def createUserInterface( self ):
    '''
    '''
    self.__layout = self.__parent.createUserInterface()
    
    loadDataButton = qt.QPushButton('Download test MRA and MRI data')
    self.__layout.addRow(loadDataButton)
    loadDataButton.connect('clicked()', self.loadData)

    inputMriLabel = qt.QLabel( 'LGE-MRI Image:' )
    self.__inputMriSelector = slicer.qMRMLNodeComboBox()
    self.__inputMriSelector.toolTip = "Select the LGE-MRI image to be analyzed"
    self.__inputMriSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    self.__inputMriSelector.setMRMLScene(slicer.mrmlScene)
    self.__inputMriSelector.addEnabled = False
    self.__inputMriSelector.renameEnabled = True
    self.__layout.addRow( inputMriLabel, self.__inputMriSelector )
    
    inputMraLabel = qt.QLabel( 'MRA Image:' )
    self.__inputMraSelector = slicer.qMRMLNodeComboBox()
    self.__inputMraSelector.toolTip = "Select the MRA image to be analyzed"
    self.__inputMraSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    self.__inputMraSelector.setMRMLScene(slicer.mrmlScene)
    self.__inputMraSelector.addEnabled = False
    self.__inputMraSelector.renameEnabled = True
    self.__layout.addRow( inputMraLabel, self.__inputMraSelector )
    
    #self.updateWidgetFromParameters(self.parameterNode())

  def loadData(self):
    vl = slicer.modules.volumes.logic()
    # Link to image on slicer wiki is found through 'copy link location' 
    vol_mri = vl.AddArchetypeVolume('http://www.na-mic.org/Wiki/images/8/8d/LAWorkflow_LGE.nrrd', 'LGE-MRI', 0)
    if vol_mri != None:
      Helper.SetVolume(vol_mri.GetID())
    vol_mra = vl.AddArchetypeVolume('http://www.na-mic.org/Wiki/images/3/3d/LAWorkflow_MRA.nrrd', 'MRA', 0)
    if vol_mra != None:
        Helper.SetVolume(vol_mra.GetID())
    
  def validate( self, desiredBranchId ):    
    self.__parent.validate( desiredBranchId )
    inputMri = self.__inputMriSelector.currentNode()
    inputMra = self.__inputMraSelector.currentNode()
    
    if inputMri != None and inputMra != None:
      inputMriID = inputMri.GetID()
      pNode = self.parameterNode()
      pNode.SetParameter('inputMriID', inputMriID)
      inputMraID = inputMra.GetID()
      pNode = self.parameterNode()
      pNode.SetParameter('inputMraID', inputMraID)
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Select Data','Please select both the images to proceed.')    

  def onEntry(self, comingFrom, transitionType):
    super(LASegmentationWorkflowSelectDataStep, self).onEntry(comingFrom, transitionType)
    #self.updateWidgetFromParameters(self.parameterNode())
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

    qt.QTimer.singleShot(0, self.killButton)

  def onExit(self, goingTo, transitionType):
    pNode = self.parameterNode()    
    #if goingTo.id() != 'LAEndoSegmentation':
    #  return
      
    super(LASegmentationWorkflowSelectDataStep, self).onExit(goingTo, transitionType)

  #def updateWidgetFromParameters(self, parameterNode):    
  #  inputVolumeID = parameterNode.GetParameter('inputVolumeID')
  #  if inputVolumeID != None:
  #    self.__inputVolumeSelector.setCurrentNode(Helper.getNodeByID(inputVolumeID))
    
    
    
