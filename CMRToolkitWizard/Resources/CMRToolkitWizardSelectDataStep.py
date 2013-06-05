from __main__ import qt, ctk

from CMRToolkitWizardStep import *
from Helper import *

class CMRToolkitWizardSelectDataStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '1. Select input data' )
    self.setDescription( 'Select Cardiac MRI data to be analyzed. Use the test data or your own images.' )

    self.__parent = super( CMRToolkitWizardSelectDataStep, self )

  def createUserInterface( self ):
    '''
    '''
    self.__layout = self.__parent.createUserInterface()
    
    loadDataButton = qt.QPushButton('Download test data')
    self.__layout.addRow(loadDataButton)
    loadDataButton.connect('clicked()', self.loadData)

    inputVolumeLabel = qt.QLabel( 'MRI Image:' )
    self.__inputVolumeSelector = slicer.qMRMLNodeComboBox()
    self.__inputVolumeSelector.toolTip = "Select the MRI image to be analyzed"
    self.__inputVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    self.__inputVolumeSelector.setMRMLScene(slicer.mrmlScene)
    self.__inputVolumeSelector.addEnabled = False
    self.__inputVolumeSelector.renameEnabled = True
    self.__layout.addRow( inputVolumeLabel, self.__inputVolumeSelector )
    
    #self.updateWidgetFromParameters(self.parameterNode())

  def loadData(self):
    vl = slicer.modules.volumes.logic()
    # Link to image on slicer wiki is found through 'copy link location' 
    vol = vl.AddArchetypeVolume('http://www.slicer.org/slicerWiki/images/b/b2/LGE-MRI.nrrd', 'LGE-MRI', 0)
    if vol != None:
      Helper.SetVolume(vol.GetID())
    
  def validate( self, desiredBranchId ):    
    self.__parent.validate( desiredBranchId )
    inputVolume = self.__inputVolumeSelector.currentNode()
    
    if inputVolume != None :
      inputID = inputVolume.GetID()    
      pNode = self.parameterNode()
      pNode.SetParameter('inputVolumeID', inputID)        
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Please select an image to proceed.')    

  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardSelectDataStep, self).onEntry(comingFrom, transitionType)
    #self.updateWidgetFromParameters(self.parameterNode())
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

  def onExit(self, goingTo, transitionType):
    pNode = self.parameterNode()    
    if goingTo.id() != 'LAEndoSegmentation':
      return
      
    super(CMRToolkitWizardSelectDataStep, self).onExit(goingTo, transitionType)

  #def updateWidgetFromParameters(self, parameterNode):    
  #  inputVolumeID = parameterNode.GetParameter('inputVolumeID')
  #  if inputVolumeID != None:
  #    self.__inputVolumeSelector.setCurrentNode(Helper.getNodeByID(inputVolumeID))
    
    
    