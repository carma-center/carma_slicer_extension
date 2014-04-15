from __main__ import vtk, qt, ctk, slicer

import LASegmentationResources
import PythonQt
#
# LASegmentationWorkflow
#

class LASegmentationWorkflow:
    def __init__(self, parent):
        parent.title = "LA Segmentation Workflow"
        parent.categories = ["Cardiac MRI Toolkit"]
        parent.dependencies = []
        parent.contributors = ["Salma Bengali (CARMA), Alan Morris (CARMA), Josh Cates (CARMA), Rob MacLeod (CARMA)"] # replace with "Firstname Lastname (Org)"
        parent.helpText = """
            This module takes the user through each step involved in segmenting LGE-MRI images by registering to gated MRA images. To use the module follow the workflow:<br><br>
            1. Select the input cardiac gated MRA image and MRI image.<br><br>2. Manually segment the left atrial (LA) endocardium region.<br><br>3. Register the MRA image to the MRI image.<br><br>4. Edit the registered segmentation as required.<br><br>More information about this module can be found at <a href=http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/>http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/</a>
            """
        parent.acknowledgementText = """
            This file was supported by...
            """ # replace with organization, grant and thanks.
        
        self.parent = parent

class LASegmentationWorkflowWidget:
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
        
        # create all 4 wizard steps
        self.selectDataStep = LASegmentationResources.LASegmentationWorkflowSelectDataStep( 'SelectData' )
        self.registrationStep = LASegmentationResources.LASegmentationWorkflowImageRegistrationStep( 'ImageRegistration' )
        self.endoSegStep = LASegmentationResources.LASegmentationWorkflowEndoSegmentationStep( 'LAEndoSegmentation' )
        
        # add the wizard steps to an array for convenience
        allSteps = []
        
        allSteps.append( self.selectDataStep )
        allSteps.append( self.registrationStep )
        allSteps.append( self.endoSegStep )        
        
        self.workflow.addTransition( self.selectDataStep, self.registrationStep )
        self.workflow.addTransition( self.registrationStep, self.endoSegStep )
        
        nNodes = slicer.mrmlScene.GetNumberOfNodesByClass('vtkMRMLScriptedModuleNode')
        
        self.parameterNode = None
        
        for n in xrange(nNodes):
            compNode = slicer.mrmlScene.GetNthNodeByClass(n, 'vtkMRMLScriptedModuleNode')
            nodeid = None
            if compNode.GetModuleName() == 'LASegmentationWorkflow':
                self.parameterNode = compNode
                print 'Found existing LASegmentationWorkflow parameter node'
                break
        
        if self.parameterNode == None:
            self.parameterNode = slicer.vtkMRMLScriptedModuleNode()
            self.parameterNode.SetModuleName('LASegmentationWorkflow')
            slicer.mrmlScene.AddNode(self.parameterNode)
        
        for s in allSteps:
            s.setParameterNode (self.parameterNode)
        
        # restore workflow step
        currentStep = self.parameterNode.GetParameter('currentStep')
        if currentStep != '':
            print 'Restoring workflow step to ', currentStep
            if currentStep == 'SelectData':
                self.workflow.setInitialStep(self.selectDataStep)
            if currentStep == 'ImageRegistration':
                self.workflow.setInitialStep(self.registrationStep)
            if currentStep == 'LAEndoSegmentation':
                self.workflow.setInitialStep(self.endoSegStep)
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
        print "LASegmentationWorkflow enter() called"
