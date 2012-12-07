from __main__ import vtk, qt, ctk, slicer

#
# CarmaRegistration
#

class CarmaRegistration:
  def __init__(self, parent):
    parent.title = "Carma Registration"
    parent.categories = ["CARMA"]
    parent.dependencies = []
    parent.contributors = ["Alan Morris (CARMA), Greg Gardner (CARMA), Salma Bengali (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This module provides registration presets for a given set of cases related to the CARMA Afib Project. It invokes the Expert Automated Registration module with different sets of parameters tuned for different registration cases.<br><br>The steps to use this module are as follows:<br><br>
1) Load the two volumes which need to be registered, and create a new output volume for the registration result.<br><br>2) Select the imaging modality based on the types of images that are being registered.<br><br>3) Mark the check box if the images have been cropped to the area around the left atrium. This will change the type of registration.<br><br>4) Click the 'Apply Registration' button to start the registration process.<br><br>5) The Advanced Registration Parameters can be modified in order to check if the registration result can be improved.<br><br>For more detailed information on these parameters, see the online documentation at: <a href=http://www.slicer.org/slicerWiki/index.php/Modules:RegisterImages-Documentation-3.6>http://www.slicer.org/slicerWiki/index.php/Modules:RegisterImages-Documentation-3.6</a><br><br>More information about this module can be found at <a href=http://www.na-mic.org/Wiki/index.php/DBP3:Utah:SlicerModuleCardiacRegistration>http://www.na-mic.org/Wiki/index.php/DBP3:Utah:SlicerModuleCardiacRegistration</a>
    """
    parent.acknowledgementText = """
    This file was supported by...
""" # replace with organization, grant and thanks.
    
    
    self.parent = parent

#
# CarmaRegistrationWidget
#

class CarmaRegistrationWidget:
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
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "CarmaRegistration Reload"
    self.layout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)

    # Register Volume collapsible button
    collapsibleButton = ctk.ctkCollapsibleButton()
    collapsibleButton.text = "Register Volume"
    self.layout.addWidget(collapsibleButton)

    # Layout within the dummy collapsible button
    formLayout = qt.QFormLayout(collapsibleButton)

    # The image volume selectors
    self.fixedFrame = qt.QFrame(collapsibleButton)
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

    self.movingFrame = qt.QFrame(collapsibleButton)
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

    self.outputFrame = qt.QFrame(collapsibleButton)
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
    self.presetFrame = qt.QFrame(collapsibleButton)
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
    self.cropCheckBox = qt.QCheckBox("Image volumes are cropped to area around LA", collapsibleButton)
    self.cropCheckBox.toolTip = "For best registration, indicate whether the volumes have been cropped to the LA ROI"
    formLayout.addWidget(self.cropCheckBox)  
    
    # Connect check box to registration parameters
    self.cropCheckBox.connect('stateChanged(int)', self.croppedImagesParameters)
  
    # Register Volume button
    registerButton = qt.QPushButton("Apply Registration")
    registerButton.toolTip = "Register Volume Two (Moving) onto Volume One (Fixed)."
    formLayout.addWidget(registerButton)
    registerButton.connect('clicked(bool)', self.onRegisterButtonClicked)
    
    
    # Advanced Registration Parameters collapsible button
    collapsibleButton2 = ctk.ctkCollapsibleButton()
    collapsibleButton2.text = "Advanced Registration Parameters"
    collapsibleButton2.collapsedHeight = 350
    collapsibleButton2.collapsed = True
    self.layout.addWidget(collapsibleButton2)
    
    
    formLayout2 = qt.QFormLayout(collapsibleButton2)    
    
    # Registration options
    self.regnFrame = qt.QFrame(collapsibleButton2)
    self.regnFrame.setLayout(qt.QHBoxLayout())
    formLayout2.addWidget(self.regnFrame)
    self.regnSelector = qt.QLabel("Select the Type of Registration: ", self.regnFrame)
    self.regnFrame.layout().addWidget(self.regnSelector)
    self.regnComboBox = qt.QComboBox()
    self.regnComboBox.addItem("Rigid")
    self.regnComboBox.addItem("Affine")
    self.regnComboBox.addItem("BSpline")
    self.regnComboBox.addItem("PipelineRigid")
    self.regnComboBox.addItem("PipelineAffine")
    self.regnComboBox.addItem("PipelineBSpline")
    self.regnComboBox.setCurrentIndex(-1)
    self.regnFrame.layout().addWidget(self.regnComboBox)
    
    # Metric options
    self.metFrame = qt.QFrame(collapsibleButton2)
    self.metFrame.setLayout(qt.QHBoxLayout())
    formLayout2.addWidget(self.metFrame)
    self.metSelector = qt.QLabel("Select the Type of Metric: ", self.metFrame)
    self.metFrame.layout().addWidget(self.metSelector)
    self.metComboBox = qt.QComboBox()
    self.metComboBox.setLayout(qt.QHBoxLayout())
    self.metComboBox.addItem("MattesMI")
    self.metComboBox.addItem("NormCorr")
    self.metComboBox.addItem("MeanSqrd")
    self.metComboBox.setCurrentIndex(-1)
    self.metFrame.layout().addWidget(self.metComboBox)
    
    # Initialization options
    self.initFrame = qt.QFrame(collapsibleButton2)
    self.initFrame.setLayout(qt.QHBoxLayout())
    formLayout2.addWidget(self.initFrame)
    self.initSelector = qt.QLabel("Select the Type of Initialization: ", self.initFrame)
    self.initFrame.layout().addWidget(self.initSelector)
    self.initComboBox = qt.QComboBox()
    self.initComboBox.setLayout(qt.QHBoxLayout())
    self.initComboBox.addItem("ImageCenters")
    self.initComboBox.addItem("CentersOfMass")
    self.initComboBox.addItem("SecondMoments")
    self.initComboBox.setCurrentIndex(-1)
    self.initFrame.layout().addWidget(self.initComboBox)
    
    # Interpolation options
    self.interpFrame = qt.QFrame(collapsibleButton2)
    self.interpFrame.setLayout(qt.QHBoxLayout())
    formLayout2.addWidget(self.interpFrame)
    self.interpSelector = qt.QLabel("Select the Type of Interpolation: ", self.interpFrame)
    self.interpFrame.layout().addWidget(self.interpSelector)
    self.interpComboBox = qt.QComboBox()
    self.interpComboBox.setLayout(qt.QHBoxLayout())
    self.interpComboBox.addItem("NearestNeighbor")
    self.interpComboBox.addItem("Linear")
    self.interpComboBox.addItem("BSpline")
    self.interpComboBox.setCurrentIndex(-1)
    self.interpFrame.layout().addWidget(self.interpComboBox)
    
    # Maximum number of iterations
    self.iterFrame = qt.QFrame(collapsibleButton2)
    self.iterFrame.setLayout(qt.QHBoxLayout())
    formLayout2.addWidget(self.iterFrame)
    self.iterEditSelector = qt.QLabel("Enter the Maximum Number of Iterations: ", self.iterFrame)
    self.iterFrame.layout().addWidget(self.iterEditSelector)
    self.iterSpinBox = qt.QSpinBox()
    self.iterSpinBox.setLayout(qt.QHBoxLayout())
    self.iterSpinBox.setMinimum(1)
    self.iterSpinBox.setMaximum(2000)
    self.iterSpinBox.setValue(0)
    self.iterFrame.layout().addWidget(self.iterSpinBox)    
    
    # Sampling Ratio value
    #self.sampFrame = qt.QFrame(collapsibleButton2)
    #self.sampFrame.setLayout(qt.QHBoxLayout())
    #formLayout2.addWidget(self.sampFrame)
    #self.sampSelector = qt.QLabel("Enter the Sampling Ratio: ", self.sampFrame)
    #self.sampFrame.layout().addWidget(self.sampSelector)
    #self.sampSpinBox = qt.QDoubleSpinBox()
    #self.sampSpinBox.setLayout(qt.QHBoxLayout())
    #self.sampSpinBox.setMinimum(0)
    #self.sampSpinBox.setMaximum(50)
    #self.sampSpinBox.setDecimals(3)
    #self.sampSpinBox.setValue(0)
    #self.sampFrame.layout().addWidget(self.sampSpinBox)        
    
    # Sampling Ratio Slider
    self.frame = qt.QFrame(collapsibleButton2)
    self.frame.setLayout(qt.QHBoxLayout())
    formLayout2.addWidget(self.frame)
    self.sampSliderSelector = qt.QLabel("Sampling Ratio (Accuracy of Registration): ", self.frame)
    #self.sampSliderSelector.tooltip = "Increasing the sampling ratio will increase registration accuracy but reduce the speed of registration"
    self.frame.layout().addWidget(self.sampSliderSelector)
    self.sampSliderFrame = ctk.ctkSliderWidget()
    self.sampSliderFrame.minimum = 0
    self.sampSliderFrame.maximum = 1
    self.sampSliderFrame.decimals = 2
    self.sampSliderFrame.singleStep = 0.01
    self.frame.layout().addWidget(self.sampSliderFrame)
       
    # Sample from fixed/moving overlap check box
    self.sampleCheckBox = qt.QCheckBox("Sample from fixed/moving image overlap", collapsibleButton2)
    self.sampleCheckBox.toolTip = "Limit metric evaluation to the fixed image region overlapped by the moving image"
    formLayout2.addWidget(self.sampleCheckBox)  
       
    # Set local var as instance attribute
    self.registerButton = registerButton


  def setParameters(self):
    method = self.presetComboBox.itemText(self.presetComboBox.currentIndex)    
    state = self.cropCheckBox.checkState()
    
    if method == 'LGE-MRI to LGE-MRI' or method == 'MRA to LGE-MRI' or method == 'CT to LGE-MRI':
      if state == 0:
        self.regnComboBox.setCurrentIndex(3)   # PipelineRigid
        self.metComboBox.setCurrentIndex(0)    # MattesMI
        self.initComboBox.setCurrentIndex(0)   # ImageCenters
        self.interpComboBox.setCurrentIndex(1) # Linear
        self.iterSpinBox.setValue(100)         
        #self.sampSpinBox.setValue(0.05)
        self.sampSliderFrame.value = 0.05  
      else:
        self.regnComboBox.setCurrentIndex(4)   # PipelineAffine
        self.metComboBox.setCurrentIndex(0)    # MattesMI
        self.initComboBox.setCurrentIndex(0)   # ImageCenters
        self.interpComboBox.setCurrentIndex(1) # Linear
        self.iterSpinBox.setValue(100)         
        #self.sampSpinBox.setValue(0.05) 
        self.sampSliderFrame.value = 0.05
    elif method == 'MRA to MRA': 
      if state == 0:
        self.regnComboBox.setCurrentIndex(3)   # PipelineRigid
        self.metComboBox.setCurrentIndex(1)    # NormCorr
        self.initComboBox.setCurrentIndex(0)   # ImageCenters
        self.interpComboBox.setCurrentIndex(1) # Linear
        self.iterSpinBox.setValue(100)
        #self.sampSpinBox.setValue(0.05)
        self.sampSliderFrame.value = 0.05
      else:
        self.regnComboBox.setCurrentIndex(4)   # PipelineAffine
        self.metComboBox.setCurrentIndex(1)    # NormCorr
        self.initComboBox.setCurrentIndex(0)   # ImageCenters
        self.interpComboBox.setCurrentIndex(1) # Linear
        self.iterSpinBox.setValue(100)
        #self.sampSpinBox.setValue(0.05)
        self.sampSliderFrame.value = 0.05
    
    
  def croppedImagesParameters(self):
    state = self.cropCheckBox.checkState()
    # CheckBox is unchecked
    if state == 0:  
      self.regnComboBox.setCurrentIndex(3) #PipelineAffine
    # CheckBox is checked
    elif state == 2:
      self.regnComboBox.setCurrentIndex(4) #PipelineRigid
      

  def onRegisterButtonClicked(self):
    print "Performing CARMA registration"
    method = self.presetComboBox.itemText(self.presetComboBox.currentIndex)
    print "method = " + method

    fixedVolume = self.fixedSelector.currentNode()
    movingVolume = self.movingSelector.currentNode()
    outputVolume = self.outputSelector.currentNode()
    
    if not (fixedVolume and movingVolume and outputVolume):
      qt.QMessageBox.critical(slicer.util.mainWindow(),
                              'Registration', 'Fixed, Moving and Output volumes are required for Registration')
      return
    
    # Get parameters to be given as inputs to Expert Automated Registration module
    param = {}
    param['fixedImage'] = fixedVolume.GetID()
    param['movingImage'] = movingVolume.GetID()
    param['resampledImage'] = outputVolume.GetID()

    # Get type of registration (Affine/Rigid/B-Spline)
    registrationType = self.regnComboBox.itemText(self.regnComboBox.currentIndex)
    
    param['registration'] = self.regnComboBox.itemText(self.regnComboBox.currentIndex)
    param['metric'] = self.metComboBox.itemText(self.metComboBox.currentIndex)
    param['initialization'] = self.initComboBox.itemText(self.initComboBox.currentIndex)
    param['interpolation'] = self.interpComboBox.itemText(self.interpComboBox.currentIndex)
    
    if registrationType == "Affine" or registrationType == "PipelineAffine" :
      param['affineMaxIterations'] = self.iterSpinBox.value
      param['affineSamplingRatio'] = self.sampSliderFrame.value
    elif registrationType == "Rigid" or registrationType == "PipelineRigid" :
      param['rigidMaxIterations'] = self.iterSpinBox.value
      param['rigidSamplingRatio'] = self.sampSliderFrame.value
    else:
      param['bsplineMaxIterations'] = self.iterSpinBox.value
      param['bsplineSamplingRatio'] = self.sampSliderFrame.value
        
    if self.sampleCheckBox.isChecked() :
      param['sampleFromOverlap'] = True
            
    slicer.cli.run( slicer.modules.expertautomatedregistration, None, param, wait_for_completion=True )
    

  def onReload(self,moduleName="CarmaRegistration"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    import imp, sys, os, slicer

    widgetName = moduleName + "Widget"

    # Reload the source code
    # - Set source file path
    # - Load the module to the global space
    filePath = eval('slicer.modules.%s.path' % moduleName.lower())
    p = os.path.dirname(filePath)
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
