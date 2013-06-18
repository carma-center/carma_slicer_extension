# Enable the 'Prefer executable CLIs' option in Slicer under Modules before running the wizard
from __main__ import qt, ctk, slicer

from CMRToolkitWizardStep import *
from Helper import *
import PythonQt

class CMRToolkitWizardEndoNVIsosurfaceStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '7. Generate Endocardium No Veins Isosurface' )
    self.setDescription( 'Select the endocardium segmentation image from which an isosurface is to be generated' )

    self.__parent = super( CMRToolkitWizardEndoNVIsosurfaceStep, self )

  def createUserInterface( self ):
    '''
    '''
    self.__layout = self.__parent.createUserInterface()
        
    #endoNVVolumeLabel = qt.QLabel( 'Endo No Veins Segmentation Image:' )
    #self.__endoNVVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__endoNVVolumeSelector.toolTip = "Select the segmentation image."
    #self.__endoNVVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__endoNVVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__endoNVVolumeSelector.addEnabled = False
    #self.__endoNVVolumeSelector.renameEnabled = True
    #self.__layout.addRow( endoNVVolumeLabel, self.__endoNVVolumeSelector )

    #self.__endoNVVolumeSelector.setCurrentNode( epiOutputVolume )
    
    #endoSegVolumeID = pNode.GetParameter('endoSegVolumeID')
    #endoSegVolume = Helper.getNodeByID(endoSegVolumeID)
    
    modelMakerButton = qt.QPushButton('Create isosurface of endocardium')
    self.__layout.addRow(modelMakerButton)
    modelMakerButton.connect('clicked()', self.startModelMaker)

  def startModelMaker(self):    
    # Load endo segmentation volume from previous step
    antrumCutVolumeID = pNode.GetParameter('antrumCutVolumeID')
    antrumCutVolume = Helper.getNodeByID(antrumCutVolumeID)
    
    # Initialize the parameters to run the Model Maker Module    
    if antrumCutVolume != None :
      param = {}
      param['InputVolume'] = antrumCutVolume.GetID()
      
      # Make a new model hierarchy node if needed
      numNodes = slicer.mrmlScene.GetNumberOfNodesByClass( "vtkMRMLModelHierarchyNode" )
      global outHierarchy 
      outHierarchy = None
      for n in xrange(numNodes):
        node = slicer.mrmlScene.GetNthNodeByClass( n, "vtkMRMLModelHierarchyNode" )
        if node.GetName() == "Endo Model":
          outHierarchy = node
          break
        
      outHierarchy = slicer.vtkMRMLModelHierarchyNode()
      outHierarchy.SetScene( slicer.mrmlScene )
      outHierarchy.SetName( "Endo Model" )
      slicer.mrmlScene.AddNode( outHierarchy )
      
      param['ModelSceneFile'] = outHierarchy
      slicer.cli.run( slicer.modules.modelmaker, None, param, wait_for_completion=True )
    else:
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Endo No Veins Isosurface', 'Endocardium segmentation image is required to run Model Maker.')
    
  def validate( self, desiredBranchId ):    
    self.__parent.validate( desiredBranchId )
    
    if outHierarchy != None:
      ModelID = outHierarchy.GetID()    
      pNode = self.parameterNode()
      pNode.SetParameter('ModelID', ModelID)        
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Please generate an isosurface to proceed.')    

  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardEndoNVIsosurfaceStep, self).onEntry(comingFrom, transitionType)
    # Make global in order to be accessible by 'startModelMaker'
    global pNode 
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

  def onExit(self, goingTo, transitionType):
    pNode = self.parameterNode()
    
    #if goingTo.id() != 'AutomaticLeftAtrialScar':
    #  return
      
    super(CMRToolkitWizardEndoNVIsosurfaceStep, self).onExit(goingTo, transitionType)

    