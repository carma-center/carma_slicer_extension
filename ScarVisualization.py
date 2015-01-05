from __main__ import vtk, qt, ctk, slicer
#
# Scar Visualization
#

class ScarVisualization:
  def __init__(self, parent):
    parent.title = "Scar Visualization"
    parent.categories = ["Cardiac MRI Toolkit"]
    parent.contributors = ["Salma Bengali (CARMA), Alan Morris (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)"]
    parent.helpText = """This module is used for scar or fibrosis visualization of the left atrium of the heart. It requires three inputs:<br>1.A cardiac LGE-MRI image<br>2.The LA endocardium segmentation<br> 3.The LA wall segmentation.<br><br> The final model can be thresholded as required to display regions of fibrosis or scar.<br><br>More information about this module can be found at <a href=http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/ScarVisualization>http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/ScarVisualization"""
    parent.acknowledgementText = """ """
    self.parent = parent

class ScarVisualizationWidget:
  def __init__(self, parent=None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()
	
  def setup(self):		
    # Input/Output collapsible button
    collapsibleButton = ctk.ctkCollapsibleButton()
    collapsibleButton.text = "IO"
    self.layout.addWidget(collapsibleButton)
		
    # Layout within the dummy collapsible button
    formLayout = qt.QFormLayout(collapsibleButton)

    # The input volume selector
    self.inputVolumeFrame = qt.QFrame(collapsibleButton)
    self.inputVolumeFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.inputVolumeFrame)
    self.inputVolumeSelector = qt.QLabel("Input Image Volume: ", self.inputVolumeFrame)
    self.inputVolumeFrame.layout().addWidget(self.inputVolumeSelector)
    self.inputVolumeSelector = slicer.qMRMLNodeComboBox(self.inputVolumeFrame)
    self.inputVolumeSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.inputVolumeSelector.addEnabled = False
    self.inputVolumeSelector.removeEnabled = True
    self.inputVolumeSelector.renameEnabled = True
    self.inputVolumeSelector.setMRMLScene( slicer.mrmlScene )
    self.inputVolumeFrame.layout().addWidget(self.inputVolumeSelector)
		
    # The endo segmentation selector
    self.endoLabelFrame = qt.QFrame(collapsibleButton)
    self.endoLabelFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.endoLabelFrame)
    self.endoLabelSelector = qt.QLabel("Endo Label Map: ", self.endoLabelFrame)
    self.endoLabelFrame.layout().addWidget(self.endoLabelSelector)
    self.endoLabelSelector = slicer.qMRMLNodeComboBox(self.endoLabelFrame)
    self.endoLabelSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.endoLabelSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 1 )
    self.endoLabelSelector.addEnabled = False
    self.endoLabelSelector.removeEnabled = True
    self.endoLabelSelector.renameEnabled = True
    self.endoLabelSelector.setMRMLScene( slicer.mrmlScene )
    self.endoLabelFrame.layout().addWidget(self.endoLabelSelector)

    # The wall segmentation selector
    self.wallLabelFrame = qt.QFrame(collapsibleButton)
    self.wallLabelFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.wallLabelFrame)
    self.wallLabelSelector = qt.QLabel("Wall Label Map: ", self.wallLabelFrame)
    self.wallLabelFrame.layout().addWidget(self.wallLabelSelector)
    self.wallLabelSelector = slicer.qMRMLNodeComboBox(self.wallLabelFrame)
    self.wallLabelSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.wallLabelSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 1 )
    self.wallLabelSelector.addEnabled = False
    self.wallLabelSelector.removeEnabled = True
    self.wallLabelSelector.renameEnabled = True
    self.wallLabelSelector.setMRMLScene( slicer.mrmlScene )
    self.wallLabelFrame.layout().addWidget(self.wallLabelSelector)
		
    # Add vertical spacer
    self.layout.addStretch(1)
		
    # Create model button
    createModelButton = qt.QPushButton("Create Endo and Wall Models")
    createModelButton.toolTip = "Create models from the LA endo and wall segmentations."
    formLayout.addWidget(createModelButton)
    createModelButton.connect('clicked(bool)', self.onCreateModelButtonClicked)
		
    # The model selector
    self.inputModelFrame = qt.QFrame(collapsibleButton)
    self.inputModelFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.inputModelFrame)
    self.inputModelSelector = qt.QLabel("Input Wall Model: ", self.inputModelFrame)
    self.inputModelFrame.layout().addWidget(self.inputModelSelector)
    self.inputModelSelector = slicer.qMRMLNodeComboBox(self.inputModelFrame)
    self.inputModelSelector.nodeTypes = ( ("vtkMRMLModelNode"), "" )
    self.inputModelSelector.addEnabled = False
    self.inputModelSelector.removeEnabled = False
    self.inputModelSelector.renameEnabled = True
    self.inputModelSelector.setMRMLScene( slicer.mrmlScene )
    self.inputModelFrame.layout().addWidget(self.inputModelSelector)
		
    # Probe volume button
    probeVolumeButton = qt.QPushButton("Probe Input Image With Wall Model")
    probeVolumeButton.toolTip = "Probe the input image volume with the wall model to map the image volume intensities to the model. The wall model will be deleted and replaced by the new output wall model."
    formLayout.addWidget(probeVolumeButton)
    probeVolumeButton.connect('clicked(bool)', self.onProbeVolumeButtonClicked)
		
    # The threshold slider
    self.thresholdFrame = qt.QFrame(collapsibleButton)
    self.thresholdFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.thresholdFrame)
    self.thresholdLabel = qt.QLabel("Threshold Range:", self.thresholdFrame)
    self.thresholdLabel.setToolTip("Set the lower and upper intensity thresholds for the model.")
    self.thresholdFrame.layout().addWidget(self.thresholdLabel)
    self.threshold = ctk.ctkRangeWidget(self.thresholdFrame)
    self.threshold.spinBoxAlignment = 0xff # Put entries on top
    self.threshold.singleStep = 0.01
    # Set a temporary minimum and maximum threshold value
    self.threshold.minimum = 0
    self.threshold.maximum = 1000
    self.threshold.singleStep = 1
    self.thresholdFrame.layout().addWidget(self.threshold)

    # Connect threshold slider to method
    self.threshold.connect( 'valuesChanged(double,double)', self.onThresholdValuesChanged )
		
    # Set local var as instance attribute
    self.createModelButton = createModelButton
    self.probeVolumeButton = probeVolumeButton
	
  def onCreateModelButtonClicked(self):
    paramEndo = {}
    paramWall = {}
    wallLabelImage = self.wallLabelSelector.currentNode()
    endoLabelImage = self.endoLabelSelector.currentNode()

    if not (endoLabelImage):
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Scar Visualization', 'Endo label image is required to generate a model.')
    else:
      paramEndo['InputVolume'] = endoLabelImage.GetID()
      # Make a new model hierarchy node if needed
      numNodes = slicer.mrmlScene.GetNumberOfNodesByClass( "vtkMRMLModelHierarchyNode" )
      global endoHierarchy 
      endoHierarchy = None
      for n in xrange(numNodes):
        node = slicer.mrmlScene.GetNthNodeByClass( n, "vtkMRMLModelHierarchyNode" )
        if node.GetName() == "Endo Model":
          endoHierarchy = node
          break
        
      endoHierarchy = slicer.vtkMRMLModelHierarchyNode()
        #  endoHierarchy = slicer.vtkMRMLModelNode()
      endoHierarchy.SetScene( slicer.mrmlScene )
      endoHierarchy.SetName( "Endo Model" )
      slicer.mrmlScene.AddNode( endoHierarchy )		
      paramEndo['ModelSceneFile'] = endoHierarchy
      paramEndo['Name'] = "Endo"
      # Pass parameters to Model Maker Module
      slicer.cli.run( slicer.modules.modelmaker, None, paramEndo, wait_for_completion=True )

    if not (wallLabelImage):
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Scar Visualization', 'Wall label image is required to generate a model.')
    else:
      paramWall['InputVolume'] = wallLabelImage.GetID()
      # Make a new model hierarchy node if needed
      numNodes = slicer.mrmlScene.GetNumberOfNodesByClass( "vtkMRMLModelHierarchyNode" )
      global wallHierarchy 
      wallHierarchy = None
      for n in xrange(numNodes):
        node = slicer.mrmlScene.GetNthNodeByClass( n, "vtkMRMLModelHierarchyNode" )
        if node.GetName() == "Wall Model":
          wallHierarchy = node
          break
        
      wallHierarchy = slicer.vtkMRMLModelHierarchyNode()
        #  wallHierarchy = slicer.vtkMRMLModelNode()
      wallHierarchy.SetScene( slicer.mrmlScene )
      wallHierarchy.SetName( "Wall Model" )
      wallHierarchy.SetModelNodeID("Wall Model")
      slicer.mrmlScene.AddNode( wallHierarchy )		
      paramWall['ModelSceneFile'] = wallHierarchy	
      paramWall['Name'] = "Wall"
      # Pass parameters to Model Maker Module
      slicer.cli.run( slicer.modules.modelmaker, None, paramWall, wait_for_completion=True )

  def onProbeVolumeButtonClicked(self): 
    # Set flag to prevent model from being thresholded without user input
    self.thresholdFlag = 0
    wallModel = self.inputModelSelector.currentNode() 
    inputVol = self.inputVolumeSelector.currentNode()
    
    self.polyData = vtk.vtkPolyData()
    if not (wallModel):
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                             'Scar Visualization', 'Wall model is required to proceed.')
    elif not (inputVol):
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                             'Scar Visualization', 'Input volume is required to proceed.')
    else:
      self.outputModel = slicer.vtkMRMLModelNode()
      self.outputModel.SetName( "Output Model" )
      slicer.mrmlScene.AddNode(self.outputModel)
      paramProbe = {}
      paramProbe['InputVolume'] = inputVol
      paramProbe['InputModel'] = wallModel
      paramProbe['OutputModel'] = self.outputModel
      # Pass parameters to Probe Volume Model
      slicer.cli.run( slicer.modules.probevolumewithmodel, None, paramProbe, wait_for_completion=True )
      # Set color of output model
      labelsColorNode = slicer.modules.colors.logic().GetColorTableNodeID(35)  # Color Warm Tint 3
      self.outputModel.GetDisplayNode().SetAndObserveColorNodeID(labelsColorNode)
      # Delete input model which is not required anymore
      slicer.mrmlScene.RemoveNode( wallModel )
      # Update threshold widget with scalar values from model
      self.polyData = self.outputModel.GetPolyData()
      self.numberOfPoints = self.polyData.GetNumberOfPoints()
      self.scalars = vtk.vtkFloatArray()
      self.scalars = self.polyData.GetPointData().GetScalars( "NRRDImage" )  
      self.range = [0] * 2
      numberOfComponents = self.scalars.GetNumberOfComponents()
      for i in xrange ( numberOfComponents ):
        self.scalars.GetRange(self.range, i)
      lowerThreshold = self.range[0]
      upperThreshold = self.range[1]
      # Set threshold range according to current model
      self.threshold.minimum = lowerThreshold
      self.threshold.maximum = upperThreshold
      self.thresholdFlag = 1
      #self.threshold.singleStep = (upperThreshold - lowerThreshold) / 1000.
  
  def onThresholdValuesChanged(self,lower,upper):
    # Check if this is the first update of the threshold range
    if self.thresholdFlag == 1:
      newScalars = vtk.vtkFloatArray()
      newScalars.SetNumberOfComponents( self.scalars.GetNumberOfComponents() )
      newScalars.SetNumberOfValues ( self.numberOfPoints )
      newScalars.SetName( "NewScalars" )
      # Set values of new scalars according to threshold
      for i in xrange( self.numberOfPoints ):
        value = self.scalars.GetValue( i )
        if value < lower:
          newScalars.SetValue( i, 0 )

        elif value > upper:
          newScalars.SetValue( i, self.range[1] )
        else:
  
          newScalars.SetValue( i, value )
      # Update the polydata with the new scalars
      self.polyData.GetPointData().AddArray( newScalars )
      self.polyData.GetPointData().RemoveArray( "NRRDImage" )
      self.polyData.Modified()
      self.outputModel.SetAndObservePolyData( self.polyData )
      self.outputModel.GetDisplayNode().SetScalarVisibility( 1 )
      self.outputModel.GetDisplayNode().SetActiveScalarName( "NewScalars" )
  
    
    
    
