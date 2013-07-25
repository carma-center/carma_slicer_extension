from __main__ import qt, ctk, slicer

from CMRToolkitWizardStep import *
from Helper import *
import PythonQt

class CMRToolkitWizardEndoSegmentationStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '2. Segmentation of endocardium' )
    self.setDescription( 'Segment the endocardium using the Editor paint tool and select the final label image to proceed.' )

    self.__parent = super( CMRToolkitWizardEndoSegmentationStep, self )

  def killButton(self):
    # Hide unneccesary button
    bl = slicer.util.findChildren(text='AutomaticLeft*')
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
    
    endoSegLabel = qt.QLabel( 'Endo Segmentation Image:' )
    self.__endoSegSelector = slicer.qMRMLNodeComboBox()
    self.__endoSegSelector.toolTip = "Choose the endo segmentation label image."
    self.__endoSegSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    self.__endoSegSelector.setMRMLScene(slicer.mrmlScene)
    self.__endoSegSelector.addEnabled = 0
    self.__layout.addRow( endoSegLabel, self.__endoSegSelector )
    
  #def startEditor(self):
  #  from Editor import EditorWidget
  #  editorWidget = EditorWidget(showVolumesFrame=True)
    
  def validate( self, desiredBranchId ):
    self.__parent.validate( desiredBranchId )
    endoSegVolume = self.__endoSegSelector.currentNode()
    
    if endoSegVolume != None :
      endoSegID = endoSegVolume.GetID()
      pNode = self.parameterNode()
      pNode.SetParameter('endoSegVolumeID', endoSegID)
      Helper.SetLabelVolume(endoSegVolume.GetID())
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Please select an endo segmentation image to proceed.')
    
  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardEndoSegmentationStep, self).onEntry(comingFrom, transitionType)
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

    qt.QTimer.singleShot(0, self.killButton)

  def onExit(self, goingTo, transitionType):    
    pNode = self.parameterNode()
    #if goingTo.id() != 'AxialDilate':
    #  return

    self.editorWidget.exit()
    
    super(CMRToolkitWizardEndoSegmentationStep, self).onExit(goingTo, transitionType) 
  
    
    
    
