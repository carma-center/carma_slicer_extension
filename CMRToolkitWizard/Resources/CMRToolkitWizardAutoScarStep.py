## Enable the 'Prefer executable CLIs' option in Slicer under Modules before running the wizard
from __main__ import qt, ctk, slicer

from CMRToolkitWizardStep import *
from Helper import *
import PythonQt

class CMRToolkitWizardAutoScarStep( CMRToolkitWizardStep ) :

  def __init__( self, stepid ):
    self.initialize( stepid )
    self.setName( '8. Run Automatic Left Atrial Scar segmentation module and generate scar isosurface' )
    self.setDescription( 'Run the automatic left atrial scar segmentation module and then run the model maker to visualize the regions of scar' )

    self.__parent = super( CMRToolkitWizardAutoScarStep, self )

  def killButton(self):
    # Hide unneccesary button
    bl = slicer.util.findChildren(text='AutomaticLeft*')
    if len(bl):
      bl[0].hide()

  def createUserInterface( self ):
    '''
    '''
    self.__layout = self.__parent.createUserInterface()
    
    ##TODO: Allow user to select the inputs or should the images be preloaded?
    
    #inputVolumeLabel = qt.QLabel( 'MRI Image:' )
    #self.__inputVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__inputVolumeSelector.toolTip = "Select the MRI image"
    #self.__inputVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__inputVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__inputVolumeSelector.addEnabled = False
    #self.__inputVolumeSelector.renameEnabled = True
    #self.__layout.addRow( inputVolumeLabel, self.__inputVolumeSelector )
    
    #wallVolumeLabel = qt.QLabel( 'Wall Segmentation Label Image:' )
    #self.__wallVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__wallVolumeSelector.toolTip = "Select the wall segmentation image"
    #self.__wallVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__wallVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__wallVolumeSelector.addEnabled = False
    #self.__wallVolumeSelector.renameEnabled = True
    #self.__layout.addRow( wallVolumeLabel, self.__wallVolumeSelector )
    
    #endoNVVolumeLabel = qt.QLabel( 'Endo No Veins Segmentation Label Image:' )
    #self.__endoNVVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__endoNVVolumeSelector.toolTip = "Select the endocardium no veins segmentation image"
    #self.__endoNVVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__endoNVVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__endoNVVolumeSelector.addEnabled = False
    #self.__endoNVVolumeSelector.renameEnabled = True
    #self.__layout.addRow( endoNVVolumeLabel, self.__endoNVVolumeSelector )
    
    #scarVolumeLabel = qt.QLabel( 'Output Scar Segmentation Image:' )
    #self.__scarVolumeSelector = slicer.qMRMLNodeComboBox()
    #self.__scarVolumeSelector.toolTip = "Create a new scar segmentation label image"
    #self.__scarVolumeSelector.nodeTypes = ['vtkMRMLScalarVolumeNode']
    #self.__scarVolumeSelector.setMRMLScene(slicer.mrmlScene)
    #self.__scarVolumeSelector.addEnabled = True
    #self.__scarVolumeSelector.renameEnabled = True
    #self.__scarVolumeSelector.baseName = "Scar Segmentation"
    #self.__layout.addRow( scarVolumeLabel, self.__epiVolumeSelector )

    autoScarButton = qt.QPushButton('Load Automatic Left Atrial Scar Module')
    self.__layout.addRow(autoScarButton)
    autoScarButton.connect('clicked()', self.startAutoScar)
    
    modelMakerButton = qt.QPushButton('Create isosurface of scar regions')
    self.__layout.addRow(modelMakerButton)
    modelMakerButton.connect('clicked()', self.startModelMaker)

  def startAutoScar(self):
    # Load epi and endo segmentation volumes from previous steps
    inputVolumeID = pNode.GetParameter('inputVolumeID')
    inputVolume = Helper.getNodeByID(inputVolumeID)
    
    wallSegVolumeID = pNode.GetParameter('wallSegVolumeID')
    wallSegVolume = Helper.getNodeByID(wallSegVolumeID)
    
    antrumCutVolumeID = pNode.GetParameter('antrumCutVolumeID')
    antrumCutVolume = Helper.getNodeByID(antrumCutVolumeID)
      
    volumesLogic = slicer.modules.volumes.logic()
    global scarOutputVolume 
    scarOutputVolume = volumesLogic.CreateAndAddLabelVolume( slicer.mrmlScene, wallSegVolume, 'Scar Segmentation Output' )
    
    if inputVolume != None and wallSegVolume != None and antrumCutVolume != None and scarOutputVolume != None:
      param = {}
      param['lgefn'] = inputVolume.GetID()
      param['lawallfn'] = wallSegVolume.GetID()
      param['laendofn'] = antrumCutVolume.GetID()
      param['outputfn'] = scarOutputVolume.GetID()
      slicer.cli.run( slicer.modules.cmrtoolkitautomaticleftatrialscar, None, param, wait_for_completion=True )
    else:
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Automatic Left Atrial Scar', 
                              'MRI, wall and endocardial no veins segmentation images are required to run Automatic Left Atrial Scar Module.')
                              
    self.updateParameters(self.parameterNode())
    
  def startModelMaker(self):    
    if scarOutputVolume != None:
      # Initialize the input parameters to run the Model Maker Module
      param = {}
      param['InputVolume'] = scarOutputVolume.GetID()
      
      # Make a new model hierarchy node if needed
      numNodes = slicer.mrmlScene.GetNumberOfNodesByClass( "vtkMRMLModelHierarchyNode" )
      outHierarchy = None
      for n in xrange(numNodes):
        node = slicer.mrmlScene.GetNthNodeByClass( n, "vtkMRMLModelHierarchyNode" )
        if node.GetName() == "Scar Model":
          outHierarchy = node
          break
        
      outHierarchy = slicer.vtkMRMLModelHierarchyNode()
      outHierarchy.SetScene( slicer.mrmlScene )
      outHierarchy.SetName( "Scar Model" )
      slicer.mrmlScene.AddNode( outHierarchy )
      
      param['ModelSceneFile'] = outHierarchy
      slicer.cli.run( slicer.modules.modelmaker, None, param, wait_for_completion=True )
    else:
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Model Maker', 
                              'Please run Automatic Left Atrial Scar Module before generating the isosurface.')

  def validate( self, desiredBranchId ):    
    self.__parent.validate( desiredBranchId )
    
    if scarOutputVolume != None:
      #scarSegID = scarOutputVolume.GetID()
      #pNode = self.parameterNode()
      #pNode.SetParameter('scarSegVolumeID', scarSegID)        
      self.__parent.validationSucceeded(desiredBranchId)
    else:
      self.__parent.validationFailed(desiredBranchId, 'Error','Automatic Scar Step did not complete successfully.')    

  def onEntry(self, comingFrom, transitionType):
    super(CMRToolkitWizardAutoScarStep, self).onEntry(comingFrom, transitionType)
    #self.updateWidgetFromParameters(self.parameterNode())
    # Make global in order to be accessible by 'startAutoScar' and 'startModelMaker'
    global pNode 
    pNode = self.parameterNode()
    pNode.SetParameter('currentStep', self.stepid)

    qt.QTimer.singleShot(0, self.killButton)
    
  def onExit(self, goingTo, transitionType):
    pNode = self.parameterNode()
    super(CMRToolkitWizardAutoScarStep, self).onExit(goingTo, transitionType)

  def updateParameters(self, parameterNode):
    # Set wall segmentation ID
    scarSegID = scarOutputVolume.GetID()
    pNode.SetParameter('scarSegVolumeID', scarSegID)
    
    
