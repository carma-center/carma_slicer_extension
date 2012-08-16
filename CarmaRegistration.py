from __main__ import vtk, qt, ctk, slicer

#
# CarmaRegistration
#

class CarmaRegistration:
  def __init__(self, parent):
    parent.title = "Carma Registration"
    parent.categories = ["CARMA"]
    parent.dependencies = []
    parent.contributors = ["Alan Morris (CARMA), Greg Gardner (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This module provides registration presets for a given set of cases related to the CARMA Afib project.<br><br>  Underneath, the Expert Automated Registration module is invoked with different sets of parameters tuned for the different registration cases.  <br><br>For more detailed information, see the online documentation at:<br><br><a href=http://www.na-mic.org/Wiki/index.php/DBP3:Utah:SlicerModuleCardiacRegistration>http://www.na-mic.org/Wiki/index.php/DBP3:Utah:SlicerModuleCardiacRegistration</a>
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

    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "CarmaRegistration Reload"
    self.layout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)


    # Collapsible button
    collapsibleButton = ctk.ctkCollapsibleButton()
    collapsibleButton.text = "Register Volume"
    self.layout.addWidget(collapsibleButton)

    # Layout within the dummy collapsible button
    formLayout = qt.QFormLayout(collapsibleButton)

    # the volume selectors
    self.fixedFrame = qt.QFrame(collapsibleButton)
    self.fixedFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.fixedFrame)
    self.fixedSelector = qt.QLabel("Fixed Volume: ", self.fixedFrame)
    self.fixedFrame.layout().addWidget(self.fixedSelector)
    self.fixedSelector = slicer.qMRMLNodeComboBox(self.fixedFrame)
    self.fixedSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.fixedSelector.addEnabled = False
    self.fixedSelector.removeEnabled = False
    self.fixedSelector.setMRMLScene( slicer.mrmlScene )
    self.fixedFrame.layout().addWidget(self.fixedSelector)

    self.movingFrame = qt.QFrame(collapsibleButton)
    self.movingFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.movingFrame)
    self.movingSelector = qt.QLabel("Moving Volume: ", self.movingFrame)
    self.movingFrame.layout().addWidget(self.movingSelector)
    self.movingSelector = slicer.qMRMLNodeComboBox(self.movingFrame)
    self.movingSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.movingSelector.addEnabled = False
    self.movingSelector.removeEnabled = False
    self.movingSelector.setMRMLScene( slicer.mrmlScene )
    self.movingFrame.layout().addWidget(self.movingSelector)

    self.outputFrame = qt.QFrame(collapsibleButton)
    self.outputFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.outputFrame)
    self.outputSelector = qt.QLabel("Output Volume: ", self.outputFrame)
    self.outputFrame.layout().addWidget(self.outputSelector)
    self.outputSelector = slicer.qMRMLNodeComboBox(self.outputFrame)
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputFrame.layout().addWidget(self.outputSelector)

    self.presetFrame = qt.QFrame(collapsibleButton)
    self.presetFrame.setLayout(qt.QHBoxLayout())
    formLayout.addWidget(self.presetFrame)
    self.presetLabel = qt.QLabel("Registration Type: ", self.presetFrame)
    self.presetFrame.layout().addWidget(self.presetLabel)

    self.presetComboBox = qt.QComboBox()
    self.presetComboBox.addItem("LGE-MRI to LGE-MRI")
    self.presetComboBox.addItem("MRA to LGE-MRI")
    self.presetComboBox.addItem("MRA to MRA")
    self.presetComboBox.addItem("CT to LGE-MRI")
    self.presetComboBox.addItem("Acute Scar to LGE-MRI")
    self.presetFrame.layout().addWidget(self.presetComboBox)

    self.cropCheckBox = qt.QCheckBox("Images are cropped to area around LA", collapsibleButton)
    self.cropCheckBox.toolTip = "For best registration, indicated whether the volumes have been croped to the LA ROI"
    formLayout.addWidget(self.cropCheckBox)


    # HelloWorld button
    registerButton = qt.QPushButton("Register Volume")
    registerButton.toolTip = "Register Moving Volume onto Fixed."
    formLayout.addWidget(registerButton)
    registerButton.connect('clicked(bool)', self.onRegisterButtonClicked)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Set local var as instance attribute
    self.registerButton = registerButton

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
    
    param = {}
    param['fixedImage'] = fixedVolume.GetID()
    param['movingImage'] = movingVolume.GetID()
    param['resampledImage'] = outputVolume.GetID()

    param['initialization'] = "ImageCenters"

    if method == 'LGE-MRI to LGE-MRI':
      param['metric'] = "NormCorr"
    else:
      param['metric'] = "MattesMI"

    if self.cropCheckBox.checked:
      param['registration'] = "PipelineAffine"
      print "cropped = true"
    else:
      param['registration'] = "PipelineRigid"
      print "cropped = false"

    slicer.cli.run( slicer.modules.expertautomatedregistration, None, param, wait_for_completion=True )



  def onReload(self,moduleName="CarmaRegistration"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    import imp, sys, os, slicer

    widgetName = moduleName + "Widget"

    # reload the source code
    # - set source file path
    # - load the module to the global space
    filePath = eval('slicer.modules.%s.path' % moduleName.lower())
    p = os.path.dirname(filePath)
    if not sys.path.__contains__(p):
      sys.path.insert(0,p)
    fp = open(filePath, "r")
    globals()[moduleName] = imp.load_module(
        moduleName, fp, filePath, ('.py', 'r', imp.PY_SOURCE))
    fp.close()

    # rebuild the widget
    # - find and hide the existing widget
    # - create a new widget in the existing parent
    parent = slicer.util.findChildren(name='%s Reload' % moduleName)[0].parent()
    for child in parent.children():
      try:
        child.hide()
      except AttributeError:
        pass
    globals()[widgetName.lower()] = eval(
        'globals()["%s"].%s(parent)' % (moduleName, widgetName))
    globals()[widgetName.lower()].setup()
