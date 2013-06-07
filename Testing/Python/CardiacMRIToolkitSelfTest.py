import os
import unittest
from __main__ import vtk, qt, ctk, slicer

class CardiacMRIToolkitSelfTest:
  def __init__(self, parent):
    parent.title = "CardiacMRIToolkitSelfTest" 
    parent.categories = ["Testing.TestCases"]
    parent.dependencies = []
    parent.contributors = ["Salma Bengali (CARMA), Alan Morris (CARMA), Brian Zenger (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)"]
    parent.helpText = """
    This module was developed as a self test for the Namic Slicer Tutorial Contest 2013    """
    parent.acknowledgementText = """
""" # replace with organization, grant and thanks.
    self.parent = parent

    # Add this test to the SelfTest module's list for discovery when the module
    # is created.  Since this module may be discovered before SelfTests itself,
    # create the list if it doesn't already exist.
    try:
      slicer.selfTests
    except AttributeError:
      slicer.selfTests = {}
    slicer.selfTests['CardiacMRIToolkitSelfTest'] = self.runTest

  def runTest(self):
    tester = CardiacMRIToolkitSelfTestCase()
    tester.runTest()
    
#
# CardiacMRIToolkitSelfTestWidget
#

class CardiacMRIToolkitSelfTestWidget:
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
    
   # Collapsible button
   testsCollapsibleButton = ctk.ctkCollapsibleButton()
   testsCollapsibleButton.text = "A collapsible button"
   self.layout.addWidget(testsCollapsibleButton)

   # Layout within the collapsible button
   formLayout = qt.QFormLayout(testsCollapsibleButton)
   
   # test buttons
   tests = ( ("Part 1: SelectMRIData",self.onSelectMRIData),("Part 2: LAEndoSegmentation", self.onLAEndoSegmentation) ) 
#,("Part 3: AxialDilate", 
#self.onAxialDilate), ("Part 3: BooleanRemove", self.onBooleanRemove), ("Part 3: LAWallCleanup", self.onLAWallCleanup), ("Part 3: PVAntrumCut",
#self.onPVAntrumCut), ("Part 3: EndoNVIsosurface", self.onEndoNVIsosurface), ("Part 3: AutomaticLeftAtrialScar", self.onAutomaticLeftAtrialScar) )

   for text,slot in tests:
      testButton = qt.QPushButton(text)
      testButton.toolTip = "Run the test."
      formLayout.addWidget(testButton)
      testButton.connect('clicked(bool)', slot)

   # Add vertical spacer
   self.layout.addStretch(1)
    
  def onSelectMRIData(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_SelectMRIData()
   
  def onLAEndoSegmentation(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_LAEndoSegmentation()

  def onAxialDilate(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_AxialDilate()

  def onBooleanRemove(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_BooleanRemove()

  def onLAWallCleanup(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_LAWallCleanup()

  def onPVAntrumCut(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_PVAntrumCut()

  def onEndoNVIsosurface(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_EndoNVIsosurface()

  def onAutomaticLeftAtrialScar(self):
   tester = CardiacMRIToolkitSelfTestCase()
   tester.setUp()
   tester.test_AutomaticLeftAtrialScar()

#
# CardiacMRIToolkitSelfTestLogic
#

class CardiacMRIToolkitSelfTestLogic:
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget
  """
  def __init__(self):
    pass

  def hasImageData(self,volumeNode):
    """This is a dummy logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      print('no volume node')
      return False
    if volumeNode.GetImageData() == None:
      print('no image data')
      return False
    return True

class CardiacMRIToolkitSelfTestCase(unittest.TestCase):
  """
  This is the test case for your scripted module.
  """

  def delayDisplay(self,message,msec=5000):
    """This utility method displays a small dialog and waits.
    This does two things: 1) it lets the event loop catch up
    to the state of the test so that rendering and widget updates
    have all taken place before the test continues and 2) it
    shows the user/developer/tester the state of the test
    so that we'll know when it breaks.
    """
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    self.delayDisplay("Closing the scene")
    layoutManager = slicer.app.layoutManager()
    layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutConventionalView)
    slicer.mrmlScene.Clear(0)

  def clickAndDrag(self,widget,button='Left',start=(10,10),end=(10,40),steps=20,modifiers=[]):
    """ Send synthetic mouse events to the specified widget (qMRMLSliceWidget or qMRMLThreeDView)
    button : "Left", "Middle", "Right", or "None"
    start, end : window coordinates for action
    steps : number of steps to move in
    modifiers : list containing zero or more of "Shift" or "Control"
    """
    style = widget.interactorStyle()
    interactor = style.GetInteractor()
    if button == 'Left':
      down = style.OnLeftButtonDown
      up = style.OnLeftButtonUp
    elif button == 'Right':
      down = style.OnRightButtonDown
      up = style.OnRightButtonUp
    elif button == 'Middle':
      down = style.OnMiddleButtonDown
      up = style.OnMiddleButtonUp
    elif button == 'None' or not button:
      down = lambda : None
      up = lambda : None
    else:
      raise Exception("Bad button - should be Left or Right, not %s" % button)
    if 'Shift' in modifiers:
      interactor.SetShiftKey(1)
    if 'Control' in modifiers:
      interactor.SetControlKey(1)
    interactor.SetEventPosition(*start)
    down()
    for step in xrange(steps):
      frac = float(step)/steps
      x = int(start[0] + frac*(end[0]-start[0]))
      y = int(start[1] + frac*(end[1]-start[1]))
      interactor.SetEventPosition(x,y)
      style.OnMouseMove()
    up()
    interactor.SetShiftKey(0)
    interactor.SetControlKey(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_SelectMRIData()
    #self.setUp()
    #self.test_LAEndoSegmentation()
    #self.setUp()
    #self.test_AxialDilate()
    #self.setUp()
    #self.test_BooleanRemove()
    #self.setUp()
    #self.test_LAWallCleanup()
    #self.setUp()
    #self.test_PVAntrumCut()
    #self.setUp()
    #self.test_EndoNVIsosurface()
    #self.setUp()
    #self.test_AutomaticLeftAtrialScar()
    
  def test_SelectMRIData(self):
    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    
    import urllib
    #downloads = (
    #    ('http://www.slicer.org/slicerWiki/images/b/b2/LGE-MRI.nrrd', 'LGE-MRI.nrrd', slicer.util.loadScene),
    #    )
    vl = slicer.modules.volumes.logic() 
    vol = vl.AddArchetypeVolume('http://www.slicer.org/slicerWiki/images/b/b2/LGE-MRI.nrrd', 'LGE-MRI', 0)
    #if vol != None:
    #  appLogic = slicer.app.applicationLogic()
    #  selectionNode = appLogic.GetSelectionNode()
    #  selectionNode.SetReferenceActiveVolumeID(vol)
    #  appLogic.PropagateVolumeSelection()
    
    #for url,name,loader in downloads:
    #  filePath = slicer.app.temporaryPath + '/' + name
    #  if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
    #    print('Requesting download %s from %s...\n' % (name, url))
    #    urllib.urlretrieve(url, filePath)
    #  if loader:
    #    print('Loading %s...\n' % (name,))
    #    loader(filePath)
    
    slicer.util.mainWindow
    self.delayDisplay('Finished with MRI Image download and loading\n')
    
    try:
      logic = CardiacMRIToolkitSelfTestLogic()
      mainWindow = slicer.util.mainWindow()
      layoutManager = slicer.app.layoutManager()
      threeDView = layoutManager.threeDWidget(0).threeDView()

      redWidget = layoutManager.sliceWidget('Red')
      redController = redWidget.sliceController()
      greenWidget = layoutManager.sliceWidget('Green')
      greenController = greenWidget.sliceController()

      self.delayDisplay('Models and Slice Model')
      mainWindow.moduleSelector().selectModule('Models')
      redWidget.sliceController().setSliceVisible(True);

      self.delayDisplay('Rotate')
      self.clickAndDrag(threeDView)

      self.delayDisplay('Zoom')
      threeDView = layoutManager.threeDWidget(0).threeDView()
      self.clickAndDrag(threeDView,button='Right')

      self.delayDisplay('Scroll Slices')
      for offset in xrange(-20,20,2):
        redController.setSliceOffsetValue(offset)

      self.delayDisplay('Test passed!')
      
    except Exception, e:
      import traceback
      traceback.print_exc()
      self.delayDisplay('Test caused exception!\n' + str(e))

  #def test_LAEndoSegmentation():
  #  self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
 #   import urllib
 #   downloads = (
 #       ('http://www.slicer.org/slicerWiki/images/f/fa/LA_Endo.nrrd', 'LA_Endo.nrrd', slicer.util.loadScene),
 #       )
    
  #  for url,name,loader in downloads:
  #    filePath = slicer.app.temporaryPath + '/' + name
  #    if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
  #      print('Requesting download %s from %s...\n' % (name, url))
  #      urllib.urlretrieve(url, filePath)
  #    if loader:
  #      print('Loading %s...\n' % (name,))
  #      loader(filePath)
  #  self.delayDisplay('Finished with LA_Endo download and loading\n')
    
  #  try:
  #    logic = CardiacMRIToolkitSelfTestLogic()
  #    mainWindow = slicer.util.mainWindow()
  #    layoutManager = slicer.app.layoutManager()
  #    threeDView = layoutManager.threeDWidget(0).threeDView()

  #    redWidget = layoutManager.sliceWidget('Red')
  #    redController = redWidget.sliceController()
  #    greenWidget = layoutManager.sliceWidget('Green')
  #    greenController = greenWidget.sliceController()

  #    self.delayDisplay('Models and Slice Model')
  #    mainWindow.moduleSelector().selectModule('Models')
  #    redWidget.sliceController().setSliceVisible(True);

  #    self.delayDisplay('Rotate')
  #    self.clickAndDrag(threeDView)

  #    self.delayDisplay('Zoom')
  #    threeDView = layoutManager.threeDWidget(0).threeDView()
  #    self.clickAndDrag(threeDView,button='Right')

  #    self.delayDisplay('Scroll Slices')
  #    for offset in xrange(-20,20,2):
  #      redController.setSliceOffsetValue(offset)

  #    self.delayDisplay('Test passed!')
      
  #  except Exception, e:
  #    import traceback
  #    traceback.print_exc()
  #    self.delayDisplay('Test caused exception!\n' + str(e))
  









