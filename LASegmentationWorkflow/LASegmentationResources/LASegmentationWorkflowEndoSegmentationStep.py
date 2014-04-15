from __main__ import qt, ctk, slicer

from LASegmentationWorkflowStep import *
from Helper import *
import PythonQt

class LASegmentationWorkflowEndoSegmentationStep( LASegmentationWorkflowStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '3. Segmentation of endocardium' )
    self.setDescription( 'Segment the endocardium using the Editor tools.' )

    self.__parent = super( LASegmentationWorkflowEndoSegmentationStep, self )

  def killButton(self):
    # Hide unneccesary button
    bl = slicer.util.findChildren(text='LAEndo*')
    if len(bl):
      bl[0].hide()

  def createUserInterface( self ):
    from Editor import EditorWidget
    
    self.__layout = self.__parent.createUserInterface()
       
    #TODO: Create label map and set it for editing by user?
    #volumesLogic = slicer.modules.volumes.logic()
    #headLabel = volumesLogic.CreateLabelVolume( slicer.mrmlScene, head, head.GetName() + '-segmentation' )
    
    #selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    #selectionNode.SetReferenceActiveVolumeID( head.GetID() )
    #selectionNode.SetReferenceActiveLabelVolumeID( headLabel.GetID() )
    #slicer.app.applicationLogic().PropagateVolumeSelection(0)
    
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
    self.editorWidget.enter()
    
  #def startEditor(self):
  #  from Editor import EditorWidget
  #  editorWidget = EditorWidget(showVolumesFrame=True)
    
  def validate( self, desiredBranchId ):
    self.__parent.validate( desiredBranchId )
    
    if endoSegVolume != None :
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Endo segmentation step did not complete successfully.')
    
  def onEntry(self, comingFrom, transitionType):
    super(LASegmentationWorkflowEndoSegmentationStep, self).onEntry(comingFrom, transitionType)
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

    qt.QTimer.singleShot(0, self.killButton)

  def onExit(self, goingTo, transitionType):    
    pNode = self.parameterNode()
    self.editorWidget.exit()
    
    super(LASegmentationWorkflowEndoSegmentationStep, self).onExit(goingTo, transitionType) 
  
    
    
    
