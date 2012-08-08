from __main__ import vtk, qt, ctk, slicer

#
# CarmaRegistration
#

class CarmaRegistration:
  def __init__(self, parent):
    parent.title = "CarmaRegistration" # TODO make this more human readable by adding spaces
    parent.categories = ["Examples"]
    parent.dependencies = []
    parent.contributors = ["Alan Morris (CARMA)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    """
    parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc. and Steve Pieper, Isomics, Inc.  and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.
    self.parent = parent

#
# qCarmaRegistrationWidget
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
    print "ok, about to perform registration"
    print "index = " + str(self.presetComboBox.currentIndex)
    print "method = " + self.presetComboBox.itemText(self.presetComboBox.currentIndex)
    if self.cropCheckBox.checked:
      print "cropped = true"
    else:
      print "cropped = false"

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
