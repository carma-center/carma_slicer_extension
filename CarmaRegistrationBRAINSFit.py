from __main__ import vtk, qt, ctk, slicer

#
# CarmaRegistrationBRAINSFit
#

class CarmaRegistrationBRAINSFit:
  def __init__(self, parent):
    parent.title = "Cardiac MRI BRAINSFit Registration"
    parent.categories = ["Cardiac MRI Toolkit"]
    parent.dependencies = []
    parent.contributors = ["Salma Bengali (CARMA), Alan Morris (CARMA), Greg Gardner (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)"]
    parent.helpText = """
    This module provides registration presets for a given set of cases related to the CARMA Afib Project. It invokes the General Registration (BRAINS) module with different sets of parameters tuned for different registration cases.<br><br>The steps to use this module are as follows:<br><br>
    1) Load the two volumes which need to be registered, and create a new output volume for the registration result.<br><br>2) Select the imaging modality based on the types of images that are being registered.<br><br>3) Mark the check box if the images have been cropped to the area around the left atrium. This will change the type of registration.<br><br>4) Click the 'Apply Registration' button to start the registration process.<br><br>5) The Advanced Registration Parameters can be modified in order to check if the registration result can be improved.<br><br>For more detailed information on these parameters, see the online documentation at: <a href=http://www.slicer.org/slicerWiki/index.php/Modules:BRAINSFit>http://www.slicer.org/slicerWiki/index.php/Modules:BRAINSFit</a><br><br>More information about this module can be found at <a href=http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/SlicerModuleCardiacRegistrationBRAINSFit>http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/SlicerModuleCardiacRegistrationBRAINSFit</a>
    """
    parent.acknowledgementText = """
    This file was supported by...
""" # replace with organization, grant and thanks.

    
    self.parent = parent
    
#
# CarmaRegistrationBRAINSFitWidget
#

class CarmaRegistrationBRAINSFitWidget:
  def __init__(self, parent = None):
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
    # Instantiate and connect widgets ...
    
    # Reload button
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "CarmaRegistrationBRAINSFit Reload"
    self.layout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)

    # Register Volume collapsible button
    registerVol = ctk.ctkCollapsibleButton()
    registerVol.text = "Register Volume"
    self.layout.addWidget(registerVol)

    # Layout within the Register Volume collapsible button
    formLayout = qt.QFormLayout(registerVol)

    # The image volume selectors
    self.fixedFrame = qt.QFrame(registerVol)
    self.fixedFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.fixedFrame)
    self.fixedSelector = qt.QLabel("Image Volume One (Fixed): ", self.fixedFrame)
    self.fixedFrame.layout().addWidget(self.fixedSelector)
    self.fixedSelector = slicer.qMRMLNodeComboBox(self.fixedFrame)
    self.fixedSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.fixedSelector.addEnabled = False
    self.fixedSelector.removeEnabled = False
    self.fixedSelector.renameEnabled = True
    self.fixedSelector.setMRMLScene( slicer.mrmlScene )
    self.fixedFrame.layout().addWidget(self.fixedSelector)

    self.movingFrame = qt.QFrame(registerVol)
    self.movingFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.movingFrame)
    self.movingSelector = qt.QLabel("Image Volume Two (Moving): ", self.movingFrame)
    self.movingFrame.layout().addWidget(self.movingSelector)
    self.movingSelector = slicer.qMRMLNodeComboBox(self.movingFrame)
    self.movingSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.movingSelector.addEnabled = False
    self.movingSelector.removeEnabled = False
    self.movingSelector.renameEnabled = True
    self.movingSelector.setMRMLScene( slicer.mrmlScene )
    self.movingFrame.layout().addWidget(self.movingSelector)

    self.outputFrame = qt.QFrame(registerVol)
    self.outputFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.outputFrame)
    self.outputSelector = qt.QLabel("Output Image Volume: ", self.outputFrame)
    self.outputFrame.layout().addWidget(self.outputSelector)
    self.outputSelector = slicer.qMRMLNodeComboBox(self.outputFrame)
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.addEnabled = True
    self.outputSelector.renameEnabled = True
    self.outputSelector.baseName = "Registered Volume"
    self.outputFrame.layout().addWidget(self.outputSelector)
    
    # Imaging Modality selection
    self.presetFrame = qt.QFrame(registerVol)
    self.presetFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.presetFrame)
    self.presetSelector = qt.QLabel("Imaging Modality: ", self.presetFrame)
    self.presetFrame.layout().addWidget(self.presetSelector)
    self.presetComboBox = qt.QComboBox()
    self.presetComboBox.addItem("LGE-MRI to LGE-MRI")
    self.presetComboBox.addItem("MRA to LGE-MRI")
    self.presetComboBox.addItem("MRA to MRA")
    self.presetComboBox.addItem("CT to LGE-MRI")
    self.presetComboBox.addItem("Acute Scar to LGE-MRI")
    self.presetComboBox.setCurrentIndex(-1)
    self.presetFrame.layout().addWidget(self.presetComboBox)    
    
    # Connect Imaging Modality widget to parameters widgets
    self.presetComboBox.connect('currentIndexChanged(QString)', self.setParameters) 
    
    # Check box
    self.cropCheckBox = qt.QCheckBox("Images are cropped to area around LA", registerVol)
    self.cropCheckBox.toolTip = "For best registration, indicate whether the volumes have been cropped to the Left Atrium ROI"
    formLayout.addWidget(self.cropCheckBox)
    
    # Connect check box to registration parameters
    self.cropCheckBox.connect('stateChanged(int)', self.croppedImagesParameters)

    # Register Volume button
    registerButton = qt.QPushButton("Apply Registration (BRAINSFit)")
    registerButton.toolTip = "Register Volume Two (Moving) onto Volume One (Fixed)."
    formLayout.addWidget(registerButton)
    registerButton.connect('clicked(bool)', self.onRegisterButtonClicked)
    
    # Initialization of Registration collapsible button
    #initOfRegn = ctk.ctkCollapsibleButton()
    #initOfRegn.text = "Initialization of Registration"
    #self.layout.addWidget(initOfRegn)

    # Layout within the Initialization of Registration collapsible button
    #formLayout2 = qt.QFormLayout(initOfRegn)
    
    # Registration Parameters collapsible button
    regnParameters = ctk.ctkCollapsibleButton()
    regnParameters.text = "Registration Parameters"
    regnParameters.collapsedHeight = 350
    regnParameters.collapsed = True
    self.layout.addWidget(regnParameters)
    
    # Layout within the Registration Parameters collapsible button
    formLayout3 = qt.QFormLayout(regnParameters)
    
    self.initOfRegnFrame = qt.QFrame(regnParameters)
    self.initOfRegnFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.initOfRegnFrame)
    self.initOfRegnSelector = qt.QLabel("Select the Type of Initialization: ", self.initOfRegnFrame)
    self.initOfRegnFrame.layout().addWidget(self.initOfRegnSelector)
    self.initOfRegnComboBox = qt.QComboBox()
    self.initOfRegnComboBox.addItem("useMomentsAlign")
    self.initOfRegnComboBox.addItem("useGeometryAlign")
    #self.initOfRegnComboBox.addItem("useCenterOfROIAlign")
    self.initOfRegnComboBox.toolTip = "GeometryAlign assumes that the center of the voxel lattice of the images represent similar structures. MomentsAlign assumes that the center of mass of the images represent similar structures."
    self.initOfRegnComboBox.setCurrentIndex(-1)
    self.initOfRegnFrame.layout().addWidget(self.initOfRegnComboBox)
    
    # Registration parameters
    self.regnFrame = qt.QFrame(regnParameters)
    self.regnFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.regnFrame)
    self.regnSelector = qt.QLabel("Select the Type of Registration: ", self.regnFrame)
    self.regnFrame.layout().addWidget(self.regnSelector)
    self.regnComboBox = qt.QComboBox()
    self.regnComboBox.addItem("Rigid")
    self.regnComboBox.addItem("Affine")
    self.regnComboBox.addItem("BSpline")
    self.regnComboBox.setCurrentIndex(-1)
    self.regnFrame.layout().addWidget(self.regnComboBox)
    
    # Maximum number of iterations
    self.iterFrame = qt.QFrame(regnParameters)
    self.iterFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.iterFrame)
    self.iterEditSelector = qt.QLabel("Enter the Maximum Number of Iterations: ", self.iterFrame)
    self.iterFrame.layout().addWidget(self.iterEditSelector)
    self.iterSpinBox = qt.QSpinBox()
    self.iterSpinBox.setLayout(qt.QHBoxLayout())
    self.iterSpinBox.setMinimum(1)
    self.iterSpinBox.setMaximum(2000)
    self.iterSpinBox.setValue(0)
    self.iterSpinBox.toolTip = "The maximum number of iterations to try before failing to converge."
    self.iterFrame.layout().addWidget(self.iterSpinBox)    

    # Number of samples
    self.sampFrame = qt.QFrame(regnParameters)
    self.sampFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.sampFrame)
    self.sampSelector = qt.QLabel("Enter the Number of Samples: ", self.sampFrame)
    self.sampFrame.layout().addWidget(self.sampSelector)
    self.sampSpinBox = qt.QSpinBox()
    self.sampSpinBox.setLayout(qt.QHBoxLayout())
    self.sampSpinBox.setMinimum(0)
    self.sampSpinBox.setMaximum(5000000)
    self.sampSpinBox.setValue(0)
    self.sampSpinBox.toolTip = "The number of voxels sampled for mutual information computation. Increase this for a slower, more careful fit."
    self.sampFrame.layout().addWidget(self.sampSpinBox)  
    
    # Minimum Step Length
    self.stepFrame = qt.QFrame(regnParameters)
    self.stepFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.stepFrame)
    self.stepSelector = qt.QLabel("Enter the Minimum Step Length: ", self.stepFrame)
    self.stepFrame.layout().addWidget(self.stepSelector)
    self.stepSpinBox = qt.QDoubleSpinBox()
    self.stepSpinBox.setLayout(qt.QHBoxLayout())
    self.stepSpinBox.setMinimum(0)
    self.stepSpinBox.setMaximum(10)
    self.stepSpinBox.setDecimals(3)
    self.stepSpinBox.setValue(0.005)
    self.stepSpinBox.toolTip = "Each step in the optimization takes steps at least this big. When none are possible, registration is complete."
    self.stepFrame.layout().addWidget(self.stepSpinBox)  
    
    # Transform Scale
    self.transFrame = qt.QFrame(regnParameters)
    self.transFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.transFrame)
    self.transSelector = qt.QLabel("Enter the Transform Scale: ", self.transFrame)
    self.transFrame.layout().addWidget(self.transSelector)
    self.transSpinBox = qt.QDoubleSpinBox()
    self.transSpinBox.setLayout(qt.QHBoxLayout())
    self.transSpinBox.setMinimum(0)
    self.transSpinBox.setMaximum(5000)
    self.transSpinBox.setValue(1000)
    self.transSpinBox.toolTip = "How much to scale up changes in position compared to unit rotational changes in radians. Decrease this to put more rotation in the search pattern."
    self.transFrame.layout().addWidget(self.transSpinBox)
    
    # Reproportion Scale
    self.repropFrame = qt.QFrame(regnParameters)
    self.repropFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.repropFrame)
    self.repropSelector = qt.QLabel("Enter the Reproportion Scale: ", self.repropFrame)
    self.repropFrame.layout().addWidget(self.repropSelector)
    self.repropSpinBox = qt.QDoubleSpinBox()
    self.repropSpinBox.setLayout(qt.QHBoxLayout())
    self.repropSpinBox.setMinimum(0)
    self.repropSpinBox.setMaximum(100)
    self.repropSpinBox.setValue(1.0)
    self.repropSpinBox.toolTip = "ScaleVersor3D 'Scale' compensation factor. Increase this to put more rescaling in a ScaleVersor3D or ScaleSkewVersor3D search pattern. 1.0 works well with a Translation Scale of 1000.0"
    self.repropFrame.layout().addWidget(self.repropSpinBox)
    
    # Skew Scale
    self.skewFrame = qt.QFrame(regnParameters)
    self.skewFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.skewFrame)
    self.skewSelector = qt.QLabel("Enter the Skew Scale: ", self.repropFrame)
    self.skewFrame.layout().addWidget(self.skewSelector)
    self.skewSpinBox = qt.QDoubleSpinBox()
    self.skewSpinBox.setLayout(qt.QHBoxLayout())
    self.skewSpinBox.setMinimum(0)
    self.skewSpinBox.setMaximum(100)
    self.skewSpinBox.setValue(1.0)
    self.skewSpinBox.toolTip = "ScaleSkewVersor3D Skew compensation factor. Increase this to put more skew in a ScaleSkewVersor3D search pattern. 1.0 works well with a translationScale of 1000.0"
    self.skewFrame.layout().addWidget(self.skewSpinBox)
    
    # Interpolation options
    self.interpFrame = qt.QFrame(regnParameters)
    self.interpFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.interpFrame)
    self.interpSelector = qt.QLabel("Select the Type of Interpolation: ", self.interpFrame)
    self.interpFrame.layout().addWidget(self.interpSelector)
    self.interpComboBox = qt.QComboBox()
    self.interpComboBox.setLayout(qt.QHBoxLayout())
    self.interpComboBox.addItem("NearestNeighbor")
    self.interpComboBox.addItem("Linear")
    self.interpComboBox.addItem("BSpline")
    self.interpComboBox.addItem("WindowedSinc")
    self.interpComboBox.addItem("RigidInPlace")
    self.interpComboBox.setCurrentIndex(-1)
    self.interpComboBox.toolTip = "Type of interpolation to be used when applying transform to moving volume."
    self.interpFrame.layout().addWidget(self.interpComboBox)
    
    # Metric Options
    self.metFrame = qt.QFrame(regnParameters)
    self.metFrame.setLayout(qt.QHBoxLayout())
    formLayout3.addWidget(self.metFrame)
    self.metSelector = qt.QLabel("Select the Type of Metric: ", self.metFrame)
    self.metFrame.layout().addWidget(self.metSelector)
    self.metComboBox = qt.QComboBox()
    self.metComboBox.setLayout(qt.QHBoxLayout())
    self.metComboBox.addItem("MMI")
    self.metComboBox.addItem("NC")
    self.metComboBox.addItem("MSE")
    self.metComboBox.addItem("MC")    
    self.metComboBox.setCurrentIndex(-1)
    self.metComboBox.toolTip = "The cost metric to be used during fitting: MMI (Mattes Mutual Information), MSE (Mean Square Error), NC (Normalized Correlation) or MC (Match Cardinality for binary images)"
    self.metFrame.layout().addWidget(self.metComboBox)


  def setParameters(self):
    method = self.presetComboBox.itemText(self.presetComboBox.currentIndex)
    state = self.cropCheckBox.checkState()

    if method == 'LGE-MRI to LGE-MRI' or method == 'MRA to LGE-MRI' or method == 'MRA to MRA':
      if state == 0:
        self.initOfRegnComboBox.setCurrentIndex(1) # useGeometryAlign
        self.regnComboBox.setCurrentIndex(0)       # Rigid
        self.iterSpinBox.setValue(200)             # Number of iterations
        self.sampSpinBox.setValue(500000)          # Number of samples
        self.stepSpinBox.setValue(0.005)           # Minimum step length
        self.transSpinBox.setValue(1000.0)         # Translation scale
        self.repropSpinBox.setValue(1.0)           # Reproportion scale
        self.skewSpinBox.setValue(1.0)             # Skew scale
        self.interpComboBox.setCurrentIndex(1)     # Linear
        self.metComboBox.setCurrentIndex(0)        # MMI
      else:
        self.initOfRegnComboBox.setCurrentIndex(1) # useGeometryAlign
        self.regnComboBox.setCurrentIndex(1)       # Affine
        self.iterSpinBox.setValue(200)             # Number of iterations
        self.sampSpinBox.setValue(500000)          # Number of samples
        self.stepSpinBox.setValue(0.005)           # Minimum step length
        self.transSpinBox.setValue(1000.0)         # Translation scale
        self.repropSpinBox.setValue(1.0)           # Reproportion scale
        self.skewSpinBox.setValue(1.0)             # Skew scale
        self.interpComboBox.setCurrentIndex(1)     # Linear
        self.metComboBox.setCurrentIndex(0)        # MMI
    elif method == 'CT to LGE-MRI':
      if state == 0:
        self.initOfRegnComboBox.setCurrentIndex(0) # useMomentsAlign
        self.regnComboBox.setCurrentIndex(0)       # Rigid
        self.iterSpinBox.setValue(200)              
        self.sampSpinBox.setValue(100000)
        self.stepSpinBox.setValue(0.005)
        self.transSpinBox.setValue(1000.0)
        self.repropSpinBox.setValue(1.0)
        self.skewSpinBox.setValue(1.0)             
        self.interpComboBox.setCurrentIndex(1)
        self.metComboBox.setCurrentIndex(0)
      else:
        self.initOfRegnComboBox.setCurrentIndex(0) # useMomentsAlign
        self.regnComboBox.setCurrentIndex(1)       # Affine
        self.iterSpinBox.setValue(200)
        self.sampSpinBox.setValue(100000)
        self.stepSpinBox.setValue(0.005)
        self.transSpinBox.setValue(1000.0)
        self.repropSpinBox.setValue(1.0)
        self.skewSpinBox.setValue(1.0)             
        self.interpComboBox.setCurrentIndex(1)
        self.metComboBox.setCurrentIndex(0)
    else:
      self.initOfRegnComboBox.setCurrentIndex(0) # useMomentsAlign
      self.regnComboBox.setCurrentIndex(0)       # Rigid
      self.iterSpinBox.setValue(200)              
      self.sampSpinBox.setValue(100000)
      self.stepSpinBox.setValue(0.005)
      self.transSpinBox.setValue(1000.0)
      self.repropSpinBox.setValue(1.0)
      self.skewSpinBox.setValue(1.0)             
      self.interpComboBox.setCurrentIndex(1)
      self.metComboBox.setCurrentIndex(0)


  def croppedImagesParameters(self):
    state = self.cropCheckBox.checkState()
    # CheckBox is unchecked
    if state == 0:  
      self.regnComboBox.setCurrentIndex(0) #Rigid
    # CheckBox is checked
    elif state == 2:
      self.regnComboBox.setCurrentIndex(1) #Affine
      

  def onRegisterButtonClicked(self):
    print "Performing CARMA registration (BRAINSFit)"
    #method = self.presetComboBox.itemText(self.presetComboBox.currentIndex)
    #print "method = " + method
    fixedVolume = self.fixedSelector.currentNode()
    movingVolume = self.movingSelector.currentNode()
    outputVolume = self.outputSelector.currentNode()
    if not (fixedVolume and movingVolume and outputVolume):
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Registration', 'Fixed, Moving and Output volumes are required for Registration')
      return
    
    param = {}
    param['fixedVolume'] = fixedVolume.GetID()
    param['movingVolume'] = movingVolume.GetID()
    param['outputVolume'] = outputVolume.GetID()
    param['transformType'] = self.regnComboBox.itemText(self.regnComboBox.currentIndex)
    param['numberOfSamples'] = self.sampSpinBox.value
    param['numberOfIterations'] = self.iterSpinBox.value
    #param['outputVolumePixelType'] = self.outputPixelComboBox.itemText(self.outputPixelComboBox.currentIndex)
    param['scaleOutputValues'] = True
    param['interpolationMode'] = self.interpComboBox.itemText(self.interpComboBox.currentIndex)
    param['minimumStepLength'] = self.stepSpinBox.value
    param['translationScale'] = self.transSpinBox.value
    param['reproportionScale'] = self.repropSpinBox.value
    param['initializeTransformMode'] = self.initOfRegnComboBox.itemText(self.initOfRegnComboBox.currentIndex)
    param['costMetric'] = self.metComboBox.itemText(self.metComboBox.currentIndex)
    param['numberOfThreads'] = 4
    
    slicer.cli.run( slicer.modules.brainsfit, None, param, wait_for_completion=True )
      
    
  def onReload(self,moduleName="CarmaRegistrationBRAINSFit"):
    
    import imp, sys, os, slicer

    widgetName = moduleName + "Widget"

    # Reload the source code
    # - Set source file path
    # - Load the module to the global space
    filePath = eval('slicer.modules.%s.path' % moduleName.lower())
    #print filePath
    p = os.path.dirname(filePath)
    #print p
    if not sys.path.__contains__(p):
      sys.path.insert(0,p)
    fp = open(filePath, "r")
    globals()[moduleName] = imp.load_module(moduleName, fp, filePath, ('.py', 'r', imp.PY_SOURCE))
    fp.close()

    # Rebuild the widget
    # - Find and hide the existing widget
    # - Create a new widget in the existing parent
    parent = slicer.util.findChildren(name='%s Reload' % moduleName)[0].parent()
    for child in parent.children():
      try:
        child.hide()
      except AttributeError:
        pass
    globals()[widgetName.lower()] = eval('globals()["%s"].%s(parent)' % (moduleName, widgetName))
    globals()[widgetName.lower()].setup()
  
