from __main__ import qt, ctk, slicer

from LASegmentationWorkflowStep import *
from Helper import *
import PythonQt

class LASegmentationWorkflowImageRegistrationStep( LASegmentationWorkflowStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '2. Image Registration' )
    self.setDescription( 'Register the MRA image to the LGE-MRI image' )

    self.__parent = super( LASegmentationWorkflowImageRegistrationStep, self )

  def killButton(self):
    # Hide unneccesary button
    bl = slicer.util.findChildren(text='LAEndo*')
    if len(bl):
      bl[0].hide()

  def createUserInterface( self ):    
    self.__layout = self.__parent.createUserInterface()
    # Register Volume button
    registerButton = qt.QPushButton("Register MRA to LGE-MRI")
    registerButton.toolTip = "Register MRA to LGE-MRI."
    self.__layout.addRow(registerButton)
    registerButton.connect('clicked(bool)', self.onRegisterButtonClicked)    

  def onRegisterButtonClicked(self):
    pNode = self.parameterNode()
    inputMriID = pNode.GetParameter('inputMriID')
    inputMri = Helper.getNodeByID(inputMriID)
    inputMraID = pNode.GetParameter('inputMraID')
    inputMra = Helper.getNodeByID(inputMraID)
    # Set parameters to be given as inputs to Expert Automated Registration module
    param = {}
    param['fixedImage'] = inputMri.GetID()
    param['movingImage'] = inputMra.GetID()
    volumesLogic = slicer.modules.volumes.logic()
    global outputVolume
    outputVolume = volumesLogic.CreateAndAddVolume( slicer.mrmlScene, inputMri, 'Registered_MRA' )
    param['resampledImage'] = outputVolume.GetID()
        
    param['registration'] = "PipelineRigid"
    param['metric'] = "MattesMI"
    param['initialization'] = "ImageCenters"
    param['interpolation'] = "Linear"
    param['rigidMaxIterations'] = 100
    param['rigidSamplingRatio'] = 0.05
    
    slicer.cli.run( slicer.modules.expertautomatedregistration, None, param, wait_for_completion=True )
    
  def validate( self, desiredBranchId ):
    self.__parent.validate( desiredBranchId )
    
    if outputVolume != None :
      outputVolumeID = outputVolume.GetID()
      pNode = self.parameterNode()
      pNode.SetParameter('outputVolumeID', outputVolumeID)
      Helper.SetLabelVolume(outputVolume.GetID())
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Registration step must be completed to proceed.')
    
  def onEntry(self, comingFrom, transitionType):
    super(LASegmentationWorkflowImageRegistrationStep, self).onEntry(comingFrom, transitionType)
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

    qt.QTimer.singleShot(0, self.killButton)

  def onExit(self, goingTo, transitionType):    
    pNode = self.parameterNode()
    
    super(LASegmentationWorkflowImageRegistrationStep, self).onExit(goingTo, transitionType) 
  
    
    
    
