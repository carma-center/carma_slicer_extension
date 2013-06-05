from __main__ import vtk, qt, ctk, slicer

import Resources
import PythonQt
#
# CMRToolkitWizard
#

class CMRToolkitWizard:
  def __init__(self, parent):
    parent.title = "CMR Toolkit Wizard"
    parent.categories = ["Cardiac MRI Toolkit"]
    parent.dependencies = []
    parent.contributors = ["Salma Bengali (CARMA), Alan Morris (CARMA), Brian Zenger (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This module takes the user through each step involved in analyzing cardiac LGE-MRI images for scar enhancement.
    """
    parent.acknowledgementText = """
    This file was supported by...
""" # replace with organization, grant and thanks.    
    
    self.parent = parent
    
#
# CMRToolkitWizardWidget
#

class CMRToolkitWizardWidget:
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
    
    if slicer.mrmlScene.GetTagByClassName( "vtkMRMLScriptedModuleNode" ) != 'ScriptedModule':
      slicer.mrmlScene.RegisterNodeClass(vtkMRMLScriptedModuleNode())
    
  def setup( self ):    
    self.workflow = ctk.ctkWorkflow()

    workflowWidget = ctk.ctkWorkflowStackedWidget()
    workflowWidget.setWorkflow( self.workflow )

    workflowWidget.buttonBoxWidget().nextButtonDefaultText = ""
    workflowWidget.buttonBoxWidget().backButtonDefaultText = ""
    
    # create all 9 wizard steps
    self.selectMRIDataStep = Resources.CMRToolkitWizardSelectDataStep( 'SelectMRIData' )
    self.endoSegStep = Resources.CMRToolkitWizardEndoSegmentationStep( 'LAEndoSegmentation' )
    self.axialDilateStep = Resources.CMRToolkitWizardAxialDilateStep( 'AxialDilate' )
    self.booleanRemoveStep = Resources.CMRToolkitWizardBooleanRemoveStep( 'BooleanRemove' )
    self.wallCleanupStep = Resources.CMRToolkitWizardWallCleanupStep( 'LAWallCleanup' )
    self.antrumCutStep = Resources.CMRToolkitWizardAntrumCutStep( 'PVAntrumCut' )
    self.endoNVIsoStep = Resources.CMRToolkitWizardEndoNVIsosurfaceStep( 'EndoNVIsosurface' )
    self.autoScarStep = Resources.CMRToolkitWizardAutoScarStep( 'AutomaticLeftAtrialScar' )
    #self.scarIsoStep = Resources.CMRToolkitWizardScarIsosurfaceStep( 'ScarIsosurface' )
    
    # add the wizard steps to an array for convenience
    allSteps = []

    allSteps.append( self.selectMRIDataStep )
    allSteps.append( self.endoSegStep )
    allSteps.append( self.axialDilateStep )
    allSteps.append( self.booleanRemoveStep )
    allSteps.append( self.wallCleanupStep )
    allSteps.append( self.antrumCutStep )
    allSteps.append( self.endoNVIsoStep )
    allSteps.append( self.autoScarStep )
    #allSteps.append( self.scarIsoStep )
    
    self.workflow.addTransition( self.selectMRIDataStep, self.endoSegStep )
    self.workflow.addTransition( self.endoSegStep, self.axialDilateStep )
    self.workflow.addTransition( self.axialDilateStep, self.booleanRemoveStep )
    self.workflow.addTransition( self.booleanRemoveStep, self.wallCleanupStep )
    self.workflow.addTransition( self.wallCleanupStep, self.antrumCutStep )
    self.workflow.addTransition( self.antrumCutStep, self.endoNVIsoStep )
    self.workflow.addTransition( self.endoNVIsoStep, self.autoScarStep )
    
    nNodes = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLScriptedModuleNode')

    self.parameterNode = None
    
    for n in xrange(nNodes):
      compNode = slicer.mrmlScene.GetNthNodeByClass(n, 'vtkMRMLScriptedModuleNode')
      nodeid = None
      if compNode.GetModuleName() == 'CMRToolkitWizard':
        self.parameterNode = compNode
        print 'Found existing CMRToolkitWizard parameter node'
        break
        
    if self.parameterNode == None:
      self.parameterNode = slicer.vtkMRMLScriptedModuleNode()
      self.parameterNode.SetModuleName('CMRToolkitWizard')
      slicer.mrmlScene.AddNode(self.parameterNode)
 
    for s in allSteps:
        s.setParameterNode (self.parameterNode)
    
     # restore workflow step
    currentStep = self.parameterNode.GetParameter('currentStep')
    if currentStep != '':
      print 'Restoring workflow step to ', currentStep
      if currentStep == 'SelectMRIData':
        self.workflow.setInitialStep(self.selectMRIDataStep)
      if currentStep == 'LAEndoSegmentation':
        self.workflow.setInitialStep(self.endoSegStep)
      if currentStep == 'AxialDilate':
        self.workflow.setInitialStep(self.axialDilateStep)
      if currentStep == 'BooleanRemove':
        self.workflow.setInitialStep(self.booleanRemoveStep)
      if currentStep == 'LAWallSegmentation':
        self.workflow.setInitialStep(self.wallCleanupStep)
      if currentStep == 'PVAntrumCut':
        self.workflow.setInitialStep(self.antrumCutStep)
      if currentStep == 'EndoNVIsosurface':
        self.workflow.setInitialStep(self.endoNVIsoStep)
      if currentStep == 'AutomaticLeftAtrialScar':
        self.workflow.setInitialStep(self.autoScarStep)
      #if currentStep == 'ScarIsosurface':
      #  self.workflow.setInitialStep(self.scarIsoStep)
    else:
      print 'currentStep in parameter node is empty!'
    
    # Start the workflow and show the widget
    self.workflow.start()
    workflowWidget.visible = True
    self.layout.addWidget( workflowWidget )

    # Enable global access to the dynamicFrames on step 2 and step 6
    #slicer.modules.emsegmentSimpleDynamicFrame = defineInputChannelsSimpleStep.dynamicFrame()
    #slicer.modules.emsegmentAdvancedDynamicFrame = definePreprocessingStep.dynamicFrame()

  def enter(self):
    print "CMRToolkitWizard enter() called"

    