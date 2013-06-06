from __main__ import qt, ctk, slicer

from CMRToolkitWizardStep import *
from Helper import *
import PythonQt

class CMRToolkitWizardWallCleanupStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '5. Edit Wall Segmentation' )
    self.setDescription( 'Set the wall segmentation image from the previous step as the merge volume in the Slicer Editor. Remove the pulmonary veins from the wall segmentation using the Slicer Editor erase tool and select the final label image to proceed.' )

    self.__parent = super( CMRToolkitWizardWallCleanupStep, self )

  def createUserInterface( self ):
    from Editor import EditorWidget
    
    self.__layout = self.__parent.createUserInterface()
    
    #endoSegmentButton = qt.QPushButton('Open Slicer Editor tool')
    #self.__layout.addRow(endoSegmentButton)
    #endoSegmentButton.connect('clicked()', self.startEditor)
    
    #TODO: Create label map and set it for editing by user?
    #volumesLogic = slicer.modules.volumes.logic()
    #headLabel = volumesLogic.CreateLabelVolume( slicer.mrmlScene, head, head.GetName() + '-segmentation' )
    
    #selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    #selectionNode.SetReferenceActiveVolumeID( head.GetID() )
    #selectionNode.SetReferenceActiveLabelVolumeID( headLabel.GetID() )
    #slicer.app.applicationLogic().PropagateVolumeSelection(0)
    
    self.updateParameters(self.parameterNode())
    
    editorFrame = qt.QFrame()
    editorFrame.setLayout(qt.QVBoxLayout())
    palette = editorFrame.palette
    bgColor = 230
    palette.setColor(qt.QPalette.Background, qt.QColor(bgColor, bgColor, bgColor))
    editorFrame.setPalette(palette)
    editorFrame.setAutoFillBackground(True);
    self.__layout.addRow(editorFrame)
    self.editorFrame = editorFrame
    global editorWidget 
    self.editorWidget = EditorWidget(parent=self.editorFrame, showVolumesFrame=True)
    self.editorWidget.setup()
    ##TODO: How to set the required merge node?
    
    wallCleanupSegLabel = qt.QLabel( 'Final Wall Segmentation Image:' )
    self.__wallCleanupSegSelector = slicer.qMRMLNodeComboBox()
    self.__wallCleanupSegSelector.toolTip = "Choose the final wall segmentation label image."
    self.__wallCleanupSegSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    self.__wallCleanupSegSelector.setMRMLScene(slicer.mrmlScene)
    self.__wallCleanupSegSelector.addEnabled = 0
    self.__layout.addRow( wallCleanupSegLabel, self.__wallCleanupSegSelector )

  def validate( self, desiredBranchId ):
    self.__parent.validate( desiredBranchId )
    wallCleanupSegVolume = self.__wallCleanupSegSelector.currentNode()
    
    if wallCleanupSegVolume != None :
      wallCleanupSegID = wallCleanupSegVolume.GetID()
      pNode = self.parameterNode()
      pNode.SetParameter('wallCleanupSegID', wallCleanupSegID)
      Helper.SetLabelVolume(wallCleanupSegVolume.GetID())
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Please select the wall segmentation image to proceed.')
    
  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardWallCleanupStep, self).onEntry(comingFrom, transitionType)
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

  ## TODO: Why does editor effect continue to the next workflow step?
  def onExit(self, goingTo, transitionType):    
    pNode = self.parameterNode()

    #if goingTo.id() != 'PVAntrumCut':
    #  return
    
    ## TODO: Why does this produce an error? 
    #self.editorWidget.exit()
    
    super(CMRToolkitWizardWallCleanupStep, self).onExit(goingTo, transitionType) 

  def updateParameters(self, parameterNode):
    global wallSegVolume 
    wallSegVolume = parameterNode.GetParameter('wallSegVolumeID')
    
    
    