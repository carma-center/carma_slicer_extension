/*
  This code performs graph cuts based on Vnet strategy, with intensity profiles as the cost functions 
*/

/*This program uses the maxflow algorithm
	described in

		"An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision."
		Yuri Boykov and Vladimir Kolmogorov.
		In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 
		September 2004

	This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
	at Siemens Corporate Research. To make it available for public use,
	it was later reimplemented by Vladimir Kolmogorov based on open publications.

	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
*/
	
/*
	For description, license, example usage see README.TXT.
*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <numeric>

#include "graph.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAddImageFilter.h"
#include "CMRToolkitLASegmentationGraphCutsCLP.h"
// Include header file for resampling image
#include "ResampleVolume.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPolygon.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataPointSampler.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryThresholdImageFunction.h"
#include "itkFloodFilledImageFunctionConditionalIterator.h"
#include "itkImageFileWriter.h"
#include "itkPluginUtilities.h"

using namespace std;

#define DATA_DMN 3 // data dimension

// defining image types
typedef itk::Image< int, DATA_DMN > ImageTypeidx;
typedef itk::Image< float,DATA_DMN > ImageType;
typedef itk::Image< unsigned char,DATA_DMN > ImageType3;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileWriter<ImageType3> WriterType;
typedef itk::AddImageFilter< ImageType, ImageType > AddImageFilterType;
typedef vector<float> vec1f;
typedef itk::Image<unsigned char, 3> LabelImageType;

#define sec_dim 25
#define stk_len 11
#define inf_vlu 1E+10
// #define DELTAS 1 // smoothness constraint
#define DELTAL 2 // lower separation constraint
#define DELTAU 12 // upper separation constraint
#define numS 2

// void printmatrix(std::vector< std::vector<int> >& matrix);
void readmatrix( std::vector< std::vector<int> >& matrix, std::ifstream& myfile, int& sz );
void connectmat( std::vector< std::vector<int> >& matrix, int sz, int** nbor_mat );
void read1Cmatrix(std::ifstream& myfile, int& num_pts);
void readCmatrix( float** &Cmatrix, std::ifstream& myfile, int& nop, vec1f& avgCOM, vec1f& realCOM);
void Mconnectmat( int** Nnbor_mat, int** nbor_mat, int& nop );
float interpolator( float* stk_index, ImageType::Pointer Image );
float dotprod( float* Mstk, float* Tstk);
void readstkmatrix( float** &stkglmatrix, std::ifstream &stkintFile, int& nop, int &numel );
float DELTASL(int f, int ALPHA, int BETA);
float DELTALB(int f, int ALPHA, int BETA);
vtkPolyData* finalmesh(float*** &Coormat, int **&Segbdrymat, float** &segmesh, std::vector< std::vector<int> >& matrix, 
           int& Mlayers, int& nop, int& sz, std::string opmeshfile, vec1f& input_origin);
void modelToLabelMap( ImageType::ConstPointer image, vtkPolyData* polydata, std::string outputImageFile );
LabelImageType::Pointer BinaryErodeFilter3D( LabelImageType::Pointer & img, unsigned int ballsize );
LabelImageType::Pointer BinaryDilateFilter3D( LabelImageType::Pointer & img, unsigned int ballsize );
LabelImageType::Pointer BinaryOpeningFilter3D( LabelImageType::Pointer & img, unsigned int ballsize );
LabelImageType::Pointer BinaryClosingFilter3D( LabelImageType::Pointer & img, unsigned int ballsize );
itk::Point<float, 3> convertStdVectorToITKPoint(const std::vector<float> & vec);

int main(int argc, char *argv[])
{
	PARSE_ARGS;
	
	std::string Coorfile, stkglfile, stkidxfile, opmeshfile1, opmeshfile2;
	//int DELTAS = atoi(argv[3]);
	int DELTAS = deltaS;
	int MDL_OPTION = MdlOption; // model selection
	// convert model option from integer to string
	std::string model_option_str;
	std::stringstream model_option;
	model_option << MDL_OPTION;
	model_option_str = model_option.str();
	vec1f RealCOM(3,0);
	
	ImageType::PointType center_of_LA;
	ImageType::IndexType center_index;
	
	vec1f MeanCOM(3,0);
	// Model's center of mass
	MeanCOM[0] = 202.088;
	MeanCOM[1] = 205.65;
	MeanCOM[2] = 57.4859;
	
	opmeshfile1 = epiMesh;
	opmeshfile2 = endoMesh;
	
	// reading input image
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inputImage);
	reader->Update();
	ImageType::ConstPointer ipimage = reader->GetOutput();
	
	// Convert point lists to ITK points and convert RAS -> LPS
	center_of_LA = convertStdVectorToITKPoint( centerOfLA );
	
	ipimage->TransformPhysicalPointToIndex( center_of_LA, center_index );
	
	RealCOM[0] = centerOfLA[0];
	RealCOM[1] = centerOfLA[1];
	RealCOM[2] = centerOfLA[2];
	
	// Get input image pixel spacing and origin
	const ImageType::SpacingType& inputspacing = ipimage->GetSpacing();
	vec1f ipspacing(3,0);
	ipspacing[0] = inputspacing[0];
	ipspacing[1] = inputspacing[1];
	ipspacing[2] = inputspacing[2];
	std::cout<< "Input image origin: " << ipimage->GetOrigin() << std::endl;
	vec1f iporigin(3,0);
	iporigin[0] = ipimage->GetOrigin()[0];
	iporigin[1] = ipimage->GetOrigin()[1];
	iporigin[2] = ipimage->GetOrigin()[2];
	
	vec1f isospacing(3,1); // isotropic voxel spacing for resampled image
	// Resampling image to isotropic pixel size
	ImageType::Pointer ResampledImage = ResampleVolumeToBe1Spacing(ipimage,ipspacing,isospacing);
	// Get Height and Width of resampled image
	//int SZ0 = ResampledImage -> GetLargestPossibleRegion().GetSize()[0];
	//int SZ1 = ResampledImage -> GetLargestPossibleRegion().GetSize()[1];
	//int SZ2 = ResampledImage -> GetLargestPossibleRegion().GetSize()[2];
	
	// Modify the centroid with respect to the resampled image
	vec1f Centroid(3,0);
	Centroid[0] = RealCOM[0] * inputspacing[0] / isospacing[0];
	Centroid[1] = RealCOM[1] * inputspacing[1] / isospacing[1];
	Centroid[2] = RealCOM[2] * inputspacing[2] / isospacing[2];
	cout << "Left Atrium center in resampled image: " << Centroid[0] << " " << Centroid[1] << " " << Centroid[2] << "\n";
	
	//file indices
	std::string model_dir = inputDataDirectory;
	
	std::string Cooridx[] = {"1_18_twin","1_17_twin","1_16_twin","1_15_twin",
		"1_14_twin","1_13_twin","1_12_twin","1_11_twin",
		"1_10_twin","1_9_twin","1_8_twin","1_7_twin","1_6_twin",
		"1_5_twin","1_4_twin","1_3_twin","1_2_twin","1_1_twin",
		"1_0_twin","0_0","0_0_twin","0_1_twin","0_2_twin",
		"0_3_twin","0_4_twin","0_5_twin","0_6_twin","0_7_twin","0_8_twin"};
	
	int layers = sizeof(Cooridx)/ sizeof( Cooridx[0] );
	std::cout << "Number of layers: "<< layers << "\n";
	//    modified layers based on inter-surface spearation constraints
	int Mlayers = layers - DELTAL;
	
	std::string textfile;
	std::stringstream triFile;
	triFile << model_dir << "/POSTENDO0" << model_option_str << "RLS.txt_0_0.tri";
	textfile = triFile.str();
	
	std::ifstream inTextFile(textfile.c_str());
	if(inTextFile.fail())
	{
		std::cerr << "Error reading connectivity file: " <<  textfile.c_str() << std::endl;
		return EXIT_FAILURE;
	}
	std::cout<<"Reading connectivity file...";
	std::vector< std::vector<int> > matrix;
	
	int sz = 0;
	readmatrix(matrix, inTextFile, sz );
	inTextFile.close();
	std::cout<<"done!"<<std::endl;
	sz--;
	
	// Initialize connectivity matrix
	int** nbor_mat = new int*[sz];
	for(int  i=0; i < sz; i++)
    nbor_mat[i] = new int[sec_dim];
	connectmat(matrix, sz, nbor_mat);
	
	// Counting the number of points of a layer
	std::stringstream PtsFile;
	PtsFile << model_dir << "/POSTENDO0" << model_option_str << "RLS.txt_0_0.pts";
	std::string PtFileStr = PtsFile.str();
	std::ifstream PtFile(PtFileStr.c_str()); // converting to C_std::string
	if(PtFile.fail())
	{
		std::cerr << "Error reading point file: " << PtFileStr.c_str() << std::endl;
		return EXIT_FAILURE;
	}
	int NumPts = 0;
	read1Cmatrix(PtFile, NumPts);
	int no_vertices = NumPts;
	std::cout << "Each layer contains: " << no_vertices << " points\n";
	
	// allocate memory for Coordinate matrix1
	float*** Coormat1;
	Coormat1 = new float**[Mlayers];
	std::cout<< "memory allocated for coordinate matrix1"<<std::endl;
	int nop = 0;
	
	// put coordinates into 3D matrix
	std::cout<< "Reading point file..." << "\n";
	for(int m = 0; m < Mlayers; m++){
		nop = 0;
		std::stringstream ptsFile;
		ptsFile << model_dir << "/POSTENDO0" << model_option_str << "RLS.txt_" << Cooridx[m] << ".pts";
		Coorfile = ptsFile.str();
		
		std::ifstream inPtFile(Coorfile.c_str()); // converting to C_std::string
		if(inPtFile.fail())
		{
			std::cerr << "Error reading point file: " << Coorfile.c_str() << std::endl;
			return EXIT_FAILURE;
		}
		Coormat1[m] = new float*[no_vertices];
		for(int i=0; i< no_vertices; i++)
			Coormat1[m][i] = new float[3];
		
		readCmatrix(Coormat1[m],inPtFile, nop, MeanCOM, Centroid);
	}
	
	// allocate memory for Coordinate matrix2
	float*** Coormat2;
	Coormat2 = new float**[Mlayers];
	std::cout<< "memory allocated for coordinate matrix2"<<std::endl;
	
	for(int m = 0; m < Mlayers; m++){
		int mDS = m + DELTAL;
		nop = 0;
		
		std::stringstream ptsFile;
		ptsFile << model_dir << "/POSTENDO0" << model_option_str << "RLS.txt_" << Cooridx[mDS] << ".pts";
		Coorfile = ptsFile.str();
		
		std::ifstream inPtFile(Coorfile.c_str()); // converting to C_std::string
		if(inPtFile.fail())
		{
			std::cerr << "Error reading point file: " << Coorfile.c_str() << std::endl;
			return EXIT_FAILURE;
		}
		Coormat2[m] = new float*[no_vertices];
		for(int i=0; i< no_vertices; i++)
			Coormat2[m][i] = new float[3];
		
		readCmatrix(Coormat2[m],inPtFile, nop, MeanCOM, Centroid);
	}
	std::cout<<"Coordinates in 3D matrix: done!"<< std::endl;
	std::cout<<"Number of points in each layer:"<< nop << std::endl;
	
	// The acquired mesh points starts from 1, so this subroutine makes it start from zero
	// by subtracting one (for programming convenience!!)
	int** Nnbor_mat = new int*[nop];
	for(int i = 0; i < nop; i++)
		Nnbor_mat[i] = new int[sec_dim];
	
	Mconnectmat(Nnbor_mat,nbor_mat, nop);
	
	// put Model stick intensities into 2D matrix
	std::stringstream OSMFile;
	OSMFile << model_dir << "/POSTWALL0" << model_option_str << "RLS_1_3_twinstk2.dat";
	stkglfile = OSMFile.str();
	
	std::ifstream stkintFile1(stkglfile.c_str()); // converting to C_std::string
	if(stkintFile1.fail())
	{
		std::cerr << "Error reading point file: " << stkglfile.c_str() << std::endl;
		return EXIT_FAILURE;
	}
	// allocate memory for Model stick intensity matrix
	int numstkgl = stk_len;
	float** Mstkglmat1;
	Mstkglmat1 = new float*[nop];
	for(int i = 0; i < nop; i++)
		Mstkglmat1[i]=new float[numstkgl];
	
	std::cout<< "memory allocated for Model1 stick intensity matrix" << std::endl;
	readstkmatrix(Mstkglmat1, stkintFile1, nop, numstkgl);
	
	std::stringstream ISMFile;
	ISMFile << model_dir << "/POSTENDO0" << model_option_str << "RLS_0_0stk2.dat";
	stkglfile = ISMFile.str();
	//stkglfile = "/Users/salmabengali/data/1/POSTENDO01RLS_0_0stk2.dat";
	
	std::ifstream stkintFile2(stkglfile.c_str()); // converting to C_std::string
	if(stkintFile2.fail())
	{
		std::cerr << "Error reading point file: " << stkglfile.c_str() << std::endl;
		return EXIT_FAILURE;
	}
	float** Mstkglmat2;
	Mstkglmat2 = new float*[nop];
	for(int i = 0; i < nop; i++)
		Mstkglmat2[i]=new float[numstkgl];
	
	std::cout<< "memory allocated for Model2 stick intensity matrix" << std::endl;
	readstkmatrix(Mstkglmat2, stkintFile2, nop, numstkgl);
	std::cout<<"Model Stick intensities in 3D matrix: done!"<<std::endl;
	
	// allocate memory for stick index matrix
	float*** stkidxmat;
	stkidxmat = new float**[layers];
	std::cout<< "memory allocated for stick index matrix" << std::endl;
	// put stick indices into 3D matrix
	for (int m = 0; m < layers; m++){
		std::stringstream MIFile;
		MIFile << model_dir << "/POSTENDO0" << model_option_str << "RLS_" << Cooridx[m] << "stkidx2.dat";;
		stkidxfile = MIFile.str();
		//stkidxfile = "/Users/salmabengali/data/1/POSTENDO01RLS_stkidx2.dat";
		//stkidxfile.insert(25 + 16,Cooridx[m]);
		
		std::ifstream stkidxFile(stkidxfile.c_str()); // converting to C_std::string
		if(stkidxFile.fail())
		{
			std::cerr << "Error reading stick/patch index file: " << stkidxfile.c_str() << std::endl;
			return EXIT_FAILURE;
		}
		int numstkidx = 3*stk_len;
		stkidxmat[m] = new float*[nop];
		for(int i = 0; i < nop; i++)
			stkidxmat[m][i]=new float[numstkidx];
		
		readstkmatrix(stkidxmat[m], stkidxFile, nop, numstkidx);
	}
	std::cout<<"Stick indices in 3D matrix: done!"<<std::endl; 
		
//    // writing Normalized cross correlation to an output file
//    std::ofstream costdata;
//    std::string costfile ="/home/sci/gveni/Documents/PhD_research/Image_segmentation/Normcostsum.dat";
//    costdata.open(costfile.c_str());
//    if ( !costdata ){
//      std::cerr << "Error: file could not be opened" << std::endl;
//      exit(1);
//    }
//    float costsum[DSsize];

    // // Reading each test image
    // for (int ii = 0; ii < 8; ii++){
    //      // if image name contains 1,2,..., use this block of code
    //      std::stringstream s;
    //      s << ii+1;
    //      std::string str = s.str();

	// allocate memory for test image stick intensity matrix
	float*** Tstkglmat1;
	Tstkglmat1 = new float**[Mlayers];
	float*** Tstkglmat2;
	Tstkglmat2 = new float**[Mlayers];
	
	for (int m = 0; m < Mlayers; m++){
		int numstkgl = stk_len;
		Tstkglmat1[m] = new float*[nop];
		Tstkglmat2[m] = new float*[nop];
		for(int i = 0; i < nop; i++){
			Tstkglmat1[m][i]=new float[numstkgl];
			Tstkglmat2[m][i]=new float[numstkgl];
		}
	}
	std::cout<< "memory allocated for Test stick intensity matrix" << std::endl;
	
	// put Test stick intensities into 3D matrix
	for (int m = 0; m < Mlayers; m++){
		for (int n = 0; n < nop; n++){
			for (int p = 0; p < stk_len; p++){
				float stk_index1[3];
				stk_index1[0] = (stkidxmat[m][n][3*p+0] - MeanCOM[0]) + Centroid[0];
				stk_index1[1] = (stkidxmat[m][n][3*p+1] - MeanCOM[1]) + Centroid[1];
				stk_index1[2] = (stkidxmat[m][n][3*p+2] - MeanCOM[2]) + Centroid[2];
				Tstkglmat1[m][n][p] = interpolator(stk_index1,ResampledImage);
				
				float stk_index2[3];
				int mDS = m + DELTAL;
				stk_index2[0] = (stkidxmat[mDS][n][3*p+0] - MeanCOM[0]) + Centroid[0];
				stk_index2[1] = (stkidxmat[mDS][n][3*p+1] - MeanCOM[1]) + Centroid[1];
				stk_index2[2] = (stkidxmat[mDS][n][3*p+2] - MeanCOM[2]) + Centroid[2];;
				Tstkglmat2[m][n][p] = interpolator(stk_index2,ResampledImage);
				
				//                 std::cout << Tstkglmat1[m][n][p] << " ";
			}
			//           std::cout << "\n";
		}
	}
	std::cout << "data entered into test stick intensity matrix" << std::endl;

	// allocate memory for corresponding cost matrix
	float** Costmat1 = new float*[Mlayers];
	float** Costmat2 = new float*[Mlayers];
	
	for(int i = 0; i < Mlayers; i++){
		Costmat1[i] = new float[nop];
		Costmat2[i] = new float[nop];
	}
	
	for (int m = 0; m < Mlayers; m++){
		for (int n = 0; n < nop; n++){
			float Mstk1[stk_len];
			float Tstk1[stk_len];
			float Mstk2[stk_len];
			float Tstk2[stk_len];
			for (int p = 0; p < stk_len; p++){
				Mstk1[p] = Mstkglmat1[n][p];
				Tstk1[p] = Tstkglmat1[m][n][p];
				Mstk2[p] = Mstkglmat2[n][p];
				Tstk2[p] = Tstkglmat2[m][n][p];
			}
			float corr1 = dotprod(Mstk1,Tstk1);
			float corr2 = dotprod(Mstk2,Tstk2);
			
			// for graph-cuts, we need cost value opposite to correlation
			float cost1,cost2;
			cost1 = -(corr1*1E4); // "without scaling, the code gives VOID results" (debugged on 7/19/2012)
			cost2 = -(corr2*1E4);
			Costmat1[m][n] = cost1;
			Costmat2[m][n] = cost2;
		}
	}
	std::cout << "costs assigned" << "\n";
	
	// allocate memory for weight matrix
	float** Wgtmat1 =  new float*[Mlayers];
	float** Wgtmat2 =  new float*[Mlayers];
	for (int i = 0; i< Mlayers; i++){
		Wgtmat1[i] = new float[nop];
		Wgtmat2[i] = new float[nop];
	}
	for (int m = 0; m < Mlayers; m++){
		for (int n = 0; n < nop; n++){
			if (m < (Mlayers-1)){
				Wgtmat1[m][n] = Costmat1[m][n] - Costmat1[m+1][n];
				Wgtmat2[m][n] = Costmat2[m][n] - Costmat2[m+1][n];
			}
			else if (m == (Mlayers-1)){
				if (n < (nop - 1)){
					Wgtmat1[m][n] = Costmat1[m][n] - (Costmat1[m][n] + Costmat1[m][n+1] - 1);
					Wgtmat2[m][n] = Costmat2[m][n] - (Costmat2[m][n] + Costmat2[m][n+1] - 1);
				}
				else{
					Wgtmat1[m][n] = Costmat1[m][n];
					Wgtmat2[m][n] = Costmat2[m][n];
				}
			}
			//            std::cout << Wgtmat1[m][n] << " ";
		}
    //          std::cout << "\n";
	}
	
	// Building a graph
	typedef Graph<float,float,float> GraphType;
	int num_nodes = Mlayers*nop;
	std::cout<< "Total number of nodes = " << num_nodes<<"\n";
	
	// calculate number of arcs from the second column of connectivity matrix
	int numedges = 0;
	for (int i = 0; i < nop; i++){
		numedges = numedges + nbor_mat[i][1];
	}
	
//  Vnet arcs
//    // num_intraarcs = T-link arcs + vertical arcs + oblique vertical arcs + base graph arcs
//    int num_intraarcs = Mlayers*nop + (Mlayers-1)*nop + (Mlayers-DELTAS)*numedges + numedges;
//    // num_interarcs = V1->V2 arcs + V2-> V1 arcs + base graph arcs
//    int num_interarcs = (Mlayers-(DELTAU-DELTAL))*nop + Mlayers*nop + numedges;
//  VCEnet arcs
    // num_intraarcs = T-link arcs + vertical arcs + oblique vertical arcs + base graph arcs
	int num_intraarcs = Mlayers*nop + (Mlayers-1)*nop + (Mlayers*(Mlayers+1)*0.5)*numedges + numedges;
	// num_interarcs = V1->V2 arcs + V2-> V1 arcs + base graph arcs
	int num_interarcs = (Mlayers-(DELTAU-DELTAL))*(DELTAU-DELTAL)*nop + Mlayers*nop + numedges;
	
	int num_arcs = numS*num_intraarcs + num_interarcs;
	std::cout<< "Total number of arcs = " << num_arcs<< "\n";
	
	GraphType *g = new GraphType(numS*num_nodes,num_arcs);
	
	// adding nodes
	std::cout << "adding nodes..." << std::endl;
	for (int gn = 0; gn < numS*num_nodes; gn++)
		g->add_node();
	
	std::cout << "nodes added. Now, adding arcs..." << std::endl;
	//adding arcs
	// arcs connecting source and sink (T-links)
	for (int m = 0; m < Mlayers; m++){
		for (int n = 0; n < nop; n++ ){
			float Tlnkcost1 = Wgtmat1[m][n];
			if(Tlnkcost1 < 0){
				g->add_tweights(m*nop + n, -Tlnkcost1, 0);
			}
			else{
				g->add_tweights(m*nop + n, 0, Tlnkcost1);
			}
			
			float Tlnkcost2 = Wgtmat2[m][n];
			if (Tlnkcost2 < 0){
				g->add_tweights(num_nodes + m*nop + n, -Tlnkcost2, 0);
			}
			else{
				g->add_tweights(num_nodes + m*nop + n, 0, Tlnkcost2);
			}
		}
	}
	std::cout << "T-links added to the graph. Now, adding N-links..." << "\n";
	
	std::cout << "arcs connecting node to the node below it"<< "\n";
	// arc from a node to the node below it
	for (int m = 0; m < (Mlayers-1); m++){
		for (int n = 0; n < nop; n++ ){
			g->add_edge(m*nop + n, ((m+1)*nop) + n, inf_vlu, 0);
			g->add_edge(num_nodes + m*nop + n, num_nodes + ((m+1)*nop) + n, inf_vlu, 0);
			//    std::cout << m * nop  + n << "<->"<< (m+1) * nop + n << " ";
		}
	}
	
	//  Vnet connections
//    std::cout<<"arcs connecting nodes in the base graph" << std::endl;
//    int m = (Mlayers-1);
//    for (int n = 0; n < nop; n++ ){
//        for (int c = 2;c < sec_dim; c++){
//            if (Nnbor_mat[n][c] != -1){
//                g->add_edge((m*nop)+Nnbor_mat[n][0],(m*nop) + Nnbor_mat[n][c], inf_vlu, 0);
//                g->add_edge(num_nodes + (m*nop)+Nnbor_mat[n][0], num_nodes + (m*nop) + Nnbor_mat[n][c], inf_vlu, 0);
//    //                    std::cout << m*nop+Nnbor_mat[n][0] << "<->" << m*nop + Nnbor_mat[n][c] << " ";
//            }
//        }
//    }

//    std::cout << "arcs connecting node to the side-nodes in the layer below it" <<std::endl;
//    for (int m = 0; m < (Mlayers-DELTAS); m++){
//        for (int n = 0; n < nop; n++ ){
//            for (int c = 2;c < sec_dim; c++){
//                if (Nnbor_mat[n][c] != -1){
//                    g->add_edge((m*nop)+Nnbor_mat[n][0],((m+DELTAS)*nop) + Nnbor_mat[n][c], inf_vlu, 0);
//                    g->add_edge(num_nodes + (m*nop)+Nnbor_mat[n][0],num_nodes + ((m+DELTAS)*nop) + Nnbor_mat[n][c], inf_vlu, 0);
//    //                            std::cout << m*nop+Nnbor_mat[n][0] << "<->" << (m+1)*nop + Nnbor_mat[n][c] << " ";
//                }
//            }
//        }
//    }

//    std::cout << "Inter-surface arcs" << "\n";
//    for (int m = 0; m < Mlayers; m++){
//        for (int n = 0; n < nop; n++ ){
//            g->add_edge(num_nodes + m*nop + n, m*nop + n, inf_vlu, 0);
//            if (m < (Mlayers-(DELTAU-DELTAL))){
//                g->add_edge(m*nop + n, num_nodes + ((m+(DELTAU-DELTAL))*nop)+n, inf_vlu, 0);
//            }
//            else if (m == (Mlayers-1))
//                g->add_edge( m*nop + n, num_nodes + m*nop + n, inf_vlu, 0);
//        }
//    }

//  VCEnet connections
    // Parameters
	int DELTAO = 4;
	int ALPHA1 = 500;
	int ALPHA2 = 5;
	int BETA = 2; // quadratic penalty
	
	//  N-links
	cout << "Oblique arcs" << "\n";
	int EdgeIntrvlIC = 10;
	vec1f Delta_ij(EdgeIntrvlIC, 0);
	for(int m = 0; m < EdgeIntrvlIC; m++){
		if(m == 0)
			Delta_ij[m] = DELTASL(m, ALPHA1, BETA);
		else
			Delta_ij[m] = DELTALB(m, ALPHA1, BETA);
	}
	for(int m = 0; m < (Mlayers-1); m++){
		for(int n = 0; n < nop; n++){
			for(int c = 2; c < sec_dim;c++){
				if(Nnbor_mat[n][c] != -1){
					for(int fij = 0; fij < EdgeIntrvlIC; fij++){
						if((m+fij) <= (Mlayers-1)){
							g->add_edge((m*nop)+Nnbor_mat[n][0],((m+fij)*nop) + Nnbor_mat[n][c],Delta_ij[fij],0);
							g->add_edge(num_nodes + (m*nop)+Nnbor_mat[n][0],num_nodes + ((m+fij)*nop) + Nnbor_mat[n][c], Delta_ij[fij], 0);
						}
					}
				}
			}
		}
	}
	cout<<"arcs connecting nodes in each base sub-graph" << endl;
	int m = (Mlayers-1);
	for (int n = 0; n < nop; n++ ){
		for (int c = 2;c < sec_dim; c++){
			if (Nnbor_mat[n][c] != -1){
				g->add_edge((m*nop)+Nnbor_mat[n][0],(m*nop) + Nnbor_mat[n][c], inf_vlu, inf_vlu);
				g->add_edge(num_nodes + (m*nop)+Nnbor_mat[n][0], num_nodes + (m*nop) + Nnbor_mat[n][c], inf_vlu, inf_vlu);
				//                    cout << m*nop+Nnbor_mat[n][0] << "<->" << m*nop + Nnbor_mat[n][c] << " ";
			}
		}
	}
	
	//  Arcs between two sub-graphs
	cout << "Inter-surface arcs" << "\n";
	int EdgeIntrvlIS = DELTAU - DELTAL;
	vec1f Delta_IS(EdgeIntrvlIS, 0);
	for(int m = 0;  m < EdgeIntrvlIS; m++){
		int q = abs((DELTAO-DELTAL) - m);
		if(m == (DELTAO-DELTAL))
			Delta_IS[m] = DELTASL(q, ALPHA2, BETA);
		else if(m != (DELTAO-DELTAL))
			Delta_IS[m] = DELTALB(q, ALPHA2, BETA);
	}
	for(int m = 0; m < Mlayers; m++){
		for(int n = 0; n < nop; n++){
			g->add_edge(num_nodes + m*nop + n, m*nop + n, inf_vlu, 0); // opposite direction arcs
			if(m < (Mlayers-1)){
				for(int fij = 1; fij < (Mlayers - m); fij++){ // fij = 1 coz no edge connection exists at fij = 0
					if((m+fij) <= (Mlayers-1))
						g->add_edge(m*nop + n, num_nodes + ((m+fij)*nop)+n, Delta_IS[fij], 0);
				}
			}else if(m == (Mlayers-1))
				g->add_edge(m*nop + n, num_nodes + m*nop + n, inf_vlu, inf_vlu);
		}
	}
	std::cout << "all arcs added" << std::endl;
	
	// Finding minimum s-t cut
	float flow = g -> maxflow();
	std::cout<< "Flow = " << flow<< "\n";
	
	// allocate memory for segmentation region matrix
	int** Segrgnmat1 =  new int*[Mlayers];
	int** Segrgnmat2 =  new int*[Mlayers];
	for (int m = 0; m< Mlayers; m++){
		Segrgnmat1[m] = new int[nop];
		Segrgnmat2[m] = new int[nop];
	}
	
	for (int m = 0; m < Mlayers; m++){
		for (int n = 0; n < nop; n++){
			if (g -> what_segment(m*nop + n) == GraphType::SOURCE)
				Segrgnmat1[m][n] = 150;
			else
				Segrgnmat1[m][n] = 0;
			
			if (g -> what_segment(num_nodes + m*nop + n) == GraphType::SOURCE)
				Segrgnmat2[m][n] = 150;
			else
				Segrgnmat2[m][n] = 0;
			//                std::cout << Segrgnmat1[m][n] << " ";
		}
    //            std::cout << "\n";
	}
	
	// allocate memory for segmentation boundary matrix
	int** Segbdrymat1 =  new int*[Mlayers];
	int** Segbdrymat2 =  new int*[Mlayers];
	for (int i = 0; i< Mlayers; i++){
		Segbdrymat1[i] = new int[nop];
		Segbdrymat2[i] = new int[nop];
	}
	
	for (int m = 0; m < Mlayers; m++){
		for (int n = 0; n < nop; n++){
			if (m == 0){
				if (Segrgnmat1[m][n] == 150)
					Segbdrymat1[m][n] = -150;
				if (Segrgnmat2[m][n] == 150)
					Segbdrymat2[m][n] = -150;
			}
			else{
				Segbdrymat1[m][n] = Segrgnmat1[m-1][n] - Segrgnmat1[m][n];
				Segbdrymat2[m][n] = Segrgnmat2[m-1][n] - Segrgnmat2[m][n];
			}
			
			//     std::cout << Segbdrymat[m][n] << " ";
		}
    //      std::cout << "\n";
	}
	
    //        //finding sum of correlation values
    //        float sumcost1 = 0;
    //        float sumcost2 = 0;

    //        for (int m = 0; m < Mlayers; m++){
    //        for (int n = 0; n < nop; n++){
    //            if (Segbdrymat1[m][n] == -150){
    //                if (Costmat1[m][n] < inf_vlu){ // To exclude divide-by-zero values
    //                    sumcost1 = sumcost1+Costmat1[m][n];
    //        //                     std::cout << sumcorr1 << " ";
    //                }
    //            }

    //           if (Segbdrymat2[m][n] == -150){
    //               if (Costmat2[m][n] < inf_vlu)
    //                   sumcost2 = sumcost2+Costmat2[m][n];
    //          }
    //        }
    //        }
    //        float sumcost = sumcost1 + sumcost2;
    //        costsum[ii] = sumcost/((float)(2*no_vertices)); // '2' coz of 2 surfaces
    //        costdata << costsum[ii] << " ";

    // build segmentation mesh
	float** segmesh1 = new float*[nop];
	float** segmesh2 = new float*[nop];
	for (int i = 0; i < nop; i++){
		segmesh1[i] = new float[3];
		segmesh2[i] = new float[3];
	}
	
	vtkPolyData *epiMeshPolyData;
	epiMeshPolyData = finalmesh(Coormat1, Segbdrymat1, segmesh1, matrix, Mlayers, nop, sz, opmeshfile1, iporigin);
	//vtkSmartPointer<vtkMRMLModelNode> epiNode = vtkSmartPointer<vtkMRMLModelNode>::New();
	//epiNode->SetAndObservePolyData( epiMeshPolyData );
	
	vtkPolyData *endoMeshPolyData; 
	endoMeshPolyData = finalmesh(Coormat2, Segbdrymat2, segmesh2, matrix, Mlayers, nop, sz, opmeshfile2, iporigin);
	
	modelToLabelMap( ipimage, endoMeshPolyData, outputEndoImage );
	modelToLabelMap( ipimage, epiMeshPolyData, outputEpiImage );		
	
	delete g; // deleting graph pointer
	return EXIT_SUCCESS;
}


void readmatrix(std::vector< std::vector<int> >& matrix, std::ifstream& myfile, int& sz)
{
  sz = 0;
  while(!myfile.eof()){
    sz = sz+1;
    std::string line;
    getline(myfile,line);
    std::stringstream lineStream(line);
    std::vector<int> numbers;
    int num;
    char strDump[5]; // for ellipsoid case
    // char strDump[1]; // for planar case

    lineStream >> strDump; // neglect the initial std::string
    lineStream >> num; // neglect
    while(lineStream >> num)
      numbers.push_back(num);

    matrix.push_back(numbers);
		//std::cerr << "Length of numbers: " << numbers.size() << std::endl;
  }
	//std::cerr << "Size of matrix: " << matrix.size() << std::endl;
}

// Design of a connectivity matrix
void connectmat(std::vector< std::vector<int> >& matrix, int sz, int** nbor_mat)
{
	
  for (int i = 0; i < sz; i++){
    for (int j= 0; j < sec_dim; j++){
      nbor_mat[i][j] = -1; // assign "-1" to the array initially
    }
  }
	
  for (int p = 0; p < sz; p++){
    int pt = p+1; // for ellipsoid case
    // int pt = p; // for planar case
    nbor_mat[p][0] = pt;
		
    int q = 1; // neighbor count
    for (int m = 0; m < sz; m++){
			
      // count # of neighbors and what are they?
      // acquire its neighbors
      if (matrix[m][0] == pt){
        q=q+1;
        nbor_mat[p][q] = matrix[m][1];
        nbor_mat[p][q+1] = matrix[m][2];
        q=q+1;
			}
			if (matrix[m][1] == pt){
				q = q+1;
				nbor_mat[p][q] = matrix[m][0];
				nbor_mat[p][q+1] = matrix[m][2];
				q=q+1;
			}
			if (matrix[m][2] == pt){
				q=q+1;
				nbor_mat[p][q] = matrix[m][0];
				nbor_mat[p][q+1] = matrix[m][1];
				q=q+1;
			}
		}
		nbor_mat[p][1] = q-1; // decrement count by 1 that adds up while counting
	}
	
  // This subroutine removes redundant neighbors and their count
  for(int i = 0; i<sz; i++){
    for (int j = 2; j<sec_dim; j++){ // leave first two columns for vertex and its neighbor count
      int curr_ele = nbor_mat[i][j];
      for (int k = j+1; k<sec_dim; k++){
        if (nbor_mat[i][k] == curr_ele && nbor_mat[i][k] != -1){
					nbor_mat[i][k] = -1;
					nbor_mat[i][1] = nbor_mat[i][1]-1;
        }
      }
		}
  }	
}

// Reading a point file to acquire number of points
void read1Cmatrix(std::ifstream& myfile, int& num_pts){
	while(!myfile.eof()){
		std::string line;
		getline(myfile,line);
		if(line.empty())
			break;
		num_pts = num_pts + 1;
	}
}

// Reading point files
void readCmatrix(float** &Cmatrix, std::ifstream& myfile, int& nop, vec1f& avgCOM, vec1f& realCOM)
{
	//std::cerr << "in readCmatrix" << std::endl;
	int v_no = 0;
	while(!myfile.eof()){
		std::string line;
		getline(myfile,line);
		if(line.empty())
			break;
		
		std::stringstream lineStream(line);
		float numbers[3];
		for(int i=0; i < 3; i++)
			lineStream >> numbers[i];
		
		Cmatrix[v_no][0] = (numbers[0]-avgCOM[0]) + realCOM[0];
		Cmatrix[v_no][1] = (numbers[1]-avgCOM[1]) + realCOM[1];
		Cmatrix[v_no][2] = (numbers[2]-avgCOM[2]) + realCOM[2];
		v_no++;
		
		nop = nop + 1;
	}
}

// stick matrix
void readstkmatrix(float** &stkglmatrix, std::ifstream &stkintFile, int& nop, int& numel){
	//std::cerr << "in readstkmatrix" << std::endl;
	for(int n = 0; n < nop; n++){
		std::string line;
		getline(stkintFile,line);
		
		std::stringstream lineStream(line);
		float *numbers = new float[numel];
		for(int i=0; i < numel; i++)
			lineStream >> numbers[i];
		
		for(int j = 0; j < numel; j++){
			stkglmatrix[n][j] = numbers[j];
			//            std::cout << stkglmatrix[n][j] << " ";
		}
		delete[] numbers;
		//        std::cout << "\n";
	}
}

// Modified connectivity matrix
void Mconnectmat(int** Nnbor_mat, int** nbor_mat, int& nop){
	for (int r = 0; r < nop; r++){
		for (int c = 0; c < sec_dim; c++){
			if (nbor_mat[r][c] != -1){
				Nnbor_mat[r][c] = nbor_mat[r][c]-1;
			}
			else
			{
				Nnbor_mat[r][c] = -1;
			}
			if (c == 1)
				Nnbor_mat[r][c] = nbor_mat[r][c]; // reverting # of neighbors to the original value
			//                std::cout<<Nnbor_mat[r][c]<< " ";
		}
		//        std::cout<<std::endl;
	}
}

float interpolator(float *stk_index, ImageType::Pointer Image){
	ImageTypeidx::IndexType VoxelIndex0;
	ImageTypeidx::IndexType VoxelIndex1;
	ImageTypeidx::IndexType VoxelIndex2;
	ImageTypeidx::IndexType VoxelIndex3;
	ImageTypeidx::IndexType VoxelIndex4;
	ImageTypeidx::IndexType VoxelIndex5;
	ImageTypeidx::IndexType VoxelIndex6;
	ImageTypeidx::IndexType VoxelIndex7;
	
	// Get size of input image
	ImageType::SizeType ImageSize = Image->GetLargestPossibleRegion().GetSize();
	
	float ptx = stk_index[0];
	float pty = stk_index[1];
	float ptz = stk_index[2];
	
	float stkglvlu;
	// Check boundary conditions
	if((ptx < 0) || (ptx > ImageSize[0]) || (pty < 0) || (pty > ImageSize[1]) || (ptz < 0) || (ptz > ImageSize[2])){
		stkglvlu = 1;
	}else{
		int ceilptx = ceil(ptx);
		int ceilpty = ceil(pty);
		int ceilptz = ceil(ptz);
		int floorptx = floor(ptx);
		int floorpty = floor(pty);
		int floorptz = floor(ptz);
		// To resolve divide-by-0 case
		if (ceilptx-floorptx == 0)
			ceilptx = ceilptx + 1;
		
		if (ceilpty-floorpty == 0)
			ceilpty = ceilpty + 1;
		
		if (ceilptz-floorptz == 0)
			ceilptz = ceilptz + 1;
		
		VoxelIndex0[0] = floorptx;
		VoxelIndex0[1] = floorpty;
		VoxelIndex0[2] = floorptz;
		ImageType::PixelType V000 = Image-> GetPixel(VoxelIndex0);
		
		VoxelIndex1[0] = floorptx;
		VoxelIndex1[1] = floorpty;
		VoxelIndex1[2] = ceilptz;
		ImageType::PixelType V001 = Image -> GetPixel(VoxelIndex1);
		
		VoxelIndex2[0] = floorptx;
		VoxelIndex2[1] = ceilpty;
		VoxelIndex2[2] = floorptz;
		ImageType::PixelType V010 = Image -> GetPixel(VoxelIndex2);
		
		VoxelIndex3[0] = floorptx;
		VoxelIndex3[1] = ceilpty;
		VoxelIndex3[2] = ceilptz;
		ImageType::PixelType V011 = Image -> GetPixel(VoxelIndex3);
		
		VoxelIndex4[0] = ceilptx;
		VoxelIndex4[1] = floorpty;
		VoxelIndex4[2] = floorptz;
		ImageType::PixelType V100 = Image -> GetPixel(VoxelIndex4);
		
		VoxelIndex5[0] = ceilptx;
		VoxelIndex5[1] = floorpty;
		VoxelIndex5[2] = ceilptz;
		ImageType::PixelType V101 = Image -> GetPixel(VoxelIndex5);
		
		VoxelIndex6[0] = ceilptx;
		VoxelIndex6[1] = ceilpty;
		VoxelIndex6[2] = floorptz;
		ImageType::PixelType V110 = Image -> GetPixel(VoxelIndex6);
		
		VoxelIndex7[0] = ceilptx;
		VoxelIndex7[1] = ceilpty;
		VoxelIndex7[2] = ceilptz;
		ImageType::PixelType V111 = Image -> GetPixel(VoxelIndex7);
		
		stkglvlu = V000 * (ceilptx-ptx) * (ceilpty-pty) * (ceilptz-ptz) +
		V100 * (ptx-floorptx) * (ceilpty-pty) * (ceilptz-ptz) +
		V010 * (ceilptx-ptx) * (pty-floorpty) * (ceilptz-ptz) +
		V001 * (ceilptx-ptx) * (ceilpty-pty) * (ptz-floorptz) +
		V101 * (ptx-floorptx) * (ceilpty-pty) * (ptz-floorptz) +
		V011 * (ceilptx-ptx) * (pty-floorpty) * (ptz-floorptz) +
		V110 * (ptx-floorptx) * (pty-floorpty) * (ceilptz-ptz) +
		V111 * (ptx-floorptx) * (pty-floorpty) * (ptz-floorptz);
	}
	return stkglvlu;
}

// compute correlation between two std::vectors
float dotprod(float *Mstk, float *Tstk){
  // Normalized correlation
	float Ncorr;
  // calculating the mean
	float Tstksum = 0;
	for (int m = 0; m < stk_len; m++){
		Tstksum = Tstksum + Tstk[m];
	}
	float Tstkmean = Tstksum/stk_len;
	
  // calculating L2 norm
	float SS_Tstk = 0;
	float varprod = 0;
	for (int m = 0; m < stk_len; m++){
		float Tvrnc = Tstk[m] - Tstkmean;
		SS_Tstk = SS_Tstk + (Tvrnc * Tvrnc);
		varprod = varprod + Mstk[m] * Tvrnc;
	}
	float Tstkl2norm = sqrt(SS_Tstk);
	if (Tstkl2norm == 0)
		Tstkl2norm = 1E5;
	
	// Normalized cross-correlation
	Ncorr = varprod/Tstkl2norm;
	return Ncorr;
}

float DELTASL(int f, int alpha, int beta){ // SL = same level
	//		double f_var = static_cast<double>(f);
	float f_var = f;
	float fa = alpha * (pow((f_var+1),beta));
	float sl = alpha * (pow(f_var,beta));
	return((fa - sl)); // cost penalty edge function
}

float DELTALB(int f, int alpha, int beta){ // LB == level below
	//		double f_var = static_cast<double>(f);
	float f_var = f;
	float fa = alpha * (pow((f_var+1),beta));
	float sl = alpha * (pow(f_var,beta));
	float fb = alpha * (pow((f_var-1),beta));
	return((fa-sl) - (sl-fb)); // cost penalty edge function
}


vtkPolyData* finalmesh(float*** &Coormat, int** &Segbdrymat, float** &segmesh, vector<vector<int> > &matrix, 
											 int &Mlayers, int &nop, int &sz, std::string opmeshfile, vec1f& input_origin){
	//std::cout << "finalmesh function" << std::endl;
	//ofstream outdata;
	//outdata.open(opmeshfile.c_str());
	//outdata.open("epi_outdata.vtk");
	/*if ( !outdata ){
	 cerr << "Error: file could not be opened" << endl;
	 return;
	 }*/
	
	//outdata << "# vtk DataFile Version 3.0\nOriginal surfacemesh\nASCII\nDATASET POLYDATA\nPOINTS " << nop << " float" << "\n";
	
	// SB: Write the points to a VTK data structure
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (int n = 0; n < nop; n++){
		for (int m = 0; m < Mlayers; m++){
			if (Segbdrymat[m][n] == -150){                
				segmesh[n][0] = -Coormat[m][n][0] - input_origin[0];
				segmesh[n][1] = -Coormat[m][n][1] - input_origin[1];
				segmesh[n][2] = +Coormat[m][n][2] + input_origin[2];
				points->InsertNextPoint( segmesh[n][0], segmesh[n][1], segmesh[n][2] );
				//outdata << segmesh[n][0] << " " << segmesh[n][1] << " " << segmesh[n][2] << "\n";
			}
		}
	}
	
	// SB: Write the polygons to a VTK data structure
	vtkSmartPointer<vtkCellArray> polygon_array = vtkSmartPointer<vtkCellArray>::New();
	//outdata << "POLYGONS " << sz << " " << 4*sz << "\n";
	for (unsigned int t = 0; t < matrix.size() - 1; t++){
		vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
		polygon->GetPointIds()->SetNumberOfIds( 3 );
		//outdata << 3;
		for (unsigned int u = 0; u < matrix[t].size(); u++){
			int val = matrix[t][u] - 1;
			polygon->GetPointIds()->SetId( u, val );
			//outdata << " " << matrix[t][u]-1;
		}
		polygon_array->InsertNextCell(polygon);
		//outdata << "\n";
	}
	//outdata.close();
	vtkSmartPointer<vtkPolyData> meshPolyData = vtkSmartPointer<vtkPolyData>::New();
	meshPolyData->SetPoints( points );
	meshPolyData->SetPolys( polygon_array );
	
	// SB: Calculate the mesh normals
	vtkSmartPointer<vtkPolyDataNormals> meshNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
	meshNormals->SetInput( meshPolyData );
	meshNormals->FlipNormalsOn();
	//meshNormals->SetFeatureAngle( 60.0 );
	//meshNormals->ComputeCellNormalsOn();
	//meshNormals->ComputePointNormalsOff();
	meshNormals->Update();
	
	// SB: Write the final mesh to a file
	vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
	writer->SetFileName(opmeshfile.c_str());
	writer->SetInputConnection(meshNormals->GetOutputPort());
	writer->Write();
	
	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData = meshNormals->GetOutput();
	return polyData;
}

void modelToLabelMap( ImageType::ConstPointer image, vtkPolyData* polyDataFile, std::string outputImageFile )
{	
	LabelImageType::Pointer label = LabelImageType::New();
  label->CopyInformation( image );
  label->SetRegions( label->GetLargestPossibleRegion() );
  label->Allocate();
  label->FillBuffer( 0 );
	
  // read the poly data
  //vtkSmartPointer<vtkPolyData> polyDataFile;
  //vtkSmartPointer<vtkPolyDataReader> pdReader;
  //vtkSmartPointer<vtkXMLPolyDataReader> pdxReader;
	
	/*pdReader = vtkSmartPointer<vtkPolyDataReader>::New();
	 pdReader->SetFileName( polyDataFile.c_str() );
	 pdReader->Update();
	 polyData = pdReader->GetOutput();*/
	
	if( polyDataFile == NULL )
	{
    std::cerr << "Failed to read model " << std::endl;
    return;
	}
	
  // LPS vs RAS
	
  vtkPoints * allPoints = polyDataFile->GetPoints();
  for( int k = 0; k < allPoints->GetNumberOfPoints(); k++ )
	{
    double* point = polyDataFile->GetPoint( k );
    point[0] = -point[0];
    point[1] = -point[1];
    allPoints->SetPoint( k, point[0], point[1], point[2] );
	}
	
  // do it
  vtkSmartPointer<vtkPolyDataPointSampler> sampler = vtkSmartPointer<vtkPolyDataPointSampler>::New();
	
  sampler->SetInput( polyDataFile );
  sampler->SetDistance( 0.5 );
  sampler->GenerateEdgePointsOn();
  sampler->GenerateInteriorPointsOn();
  sampler->GenerateVertexPointsOn();
  sampler->Update();
	
  //std::cout << polyData->GetNumberOfPoints() << std::endl;
  //std::cout << sampler->GetOutput()->GetNumberOfPoints() << std::endl;
  for( int k = 0; k < sampler->GetOutput()->GetNumberOfPoints(); k++ )
	{
    double * pt = sampler->GetOutput()->GetPoint( k );
    LabelImageType::PointType pitk;
    pitk[0] = pt[0];
    pitk[1] = pt[1];
    pitk[2] = pt[2];
    LabelImageType::IndexType idx;
    label->TransformPhysicalPointToIndex( pitk, idx );
		
    if( label->GetLargestPossibleRegion().IsInside(idx) )
		{
      label->SetPixel( idx, 255 );
		}
	}
	
  // do morphological closing
  LabelImageType::Pointer closedLabel = BinaryClosingFilter3D( label, 2);
  itk::ImageRegionIteratorWithIndex<LabelImageType> itLabel(closedLabel, closedLabel->GetLargestPossibleRegion() );
	
  // do flood fill using binary threshold image function
  typedef itk::BinaryThresholdImageFunction<LabelImageType> ImageFunctionType;
  ImageFunctionType::Pointer func = ImageFunctionType::New();
  func->SetInputImage( closedLabel );
  func->ThresholdBelow(1);
	
  LabelImageType::IndexType idx;
  LabelImageType::PointType COG;
	
  // set the centre of gravity
  // double *bounds = polyData->GetBounds();
  COG.Fill(0.0);
  for( vtkIdType k = 0; k < polyDataFile->GetNumberOfPoints(); k++ )
	{
    double *pt = polyDataFile->GetPoint( k );
    for( int m = 0; m < 3; m++ )
		{
      COG[m] += pt[m];
		}
	}
  for( int m = 0; m < 3; m++ )
	{
    COG[m] /= static_cast<float>( polyDataFile->GetNumberOfPoints() );
	}
	
  label->TransformPhysicalPointToIndex( COG, idx );
	
  itk::FloodFilledImageFunctionConditionalIterator<LabelImageType, ImageFunctionType> floodFill( closedLabel, func, idx );
  for( floodFill.GoToBegin(); !floodFill.IsAtEnd(); ++floodFill )
	{
    LabelImageType::IndexType i = floodFill.GetIndex();
    closedLabel->SetPixel( i, 255 );
	}
  LabelImageType::Pointer finalLabel = BinaryClosingFilter3D( closedLabel, 2);
  for( itLabel.GoToBegin(); !itLabel.IsAtEnd(); ++itLabel )
	{
    LabelImageType::IndexType i = itLabel.GetIndex();
    label->SetPixel( i, finalLabel->GetPixel(i) );
	}
	
	WriterType::Pointer writer = WriterType::New();
  //itk::PluginFilterWatcher watchWriter(writer, "Write Volume", CLPProcessInformation);
  writer->SetFileName( outputImageFile );
  writer->SetInput( label );
  writer->SetUseCompression(1);
  writer->Update();
}

LabelImageType::Pointer BinaryErodeFilter3D( LabelImageType::Pointer & img, unsigned int ballsize )
{
  typedef itk::BinaryBallStructuringElement<unsigned char, 3> KernalType;
  typedef itk::BinaryErodeImageFilter<LabelImageType, LabelImageType, KernalType> ErodeFilterType;
  ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
  erodeFilter->SetInput( img );

  KernalType ball;
  KernalType::SizeType ballSize;
	
  for( int k = 0; k < 3; k++ )
	{
    ballSize[k] = ballsize;
	}
	
  ball.SetRadius(ballSize);
  ball.CreateStructuringElement();
  erodeFilter->SetKernel( ball );
  erodeFilter->Update();
  return erodeFilter->GetOutput();
}

LabelImageType::Pointer BinaryDilateFilter3D( LabelImageType::Pointer & img, unsigned int ballsize )
{
  typedef itk::BinaryBallStructuringElement<unsigned char, 3> KernalType;
  typedef itk::BinaryDilateImageFilter<LabelImageType, LabelImageType, KernalType> DilateFilterType;
  DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
  dilateFilter->SetInput( img );
  KernalType ball;
  KernalType::SizeType ballSize;
	
  for( int k = 0; k < 3; k++ )
	{
    ballSize[k] = ballsize;
	}
	
  ball.SetRadius(ballSize);
  ball.CreateStructuringElement();
  dilateFilter->SetKernel( ball );
  dilateFilter->Update();
  return dilateFilter->GetOutput();
}

LabelImageType::Pointer BinaryOpeningFilter3D( LabelImageType::Pointer & img, unsigned int ballsize )
{
  LabelImageType::Pointer imgErode = BinaryErodeFilter3D( img, ballsize);

  return BinaryDilateFilter3D( imgErode, ballsize );
}

LabelImageType::Pointer BinaryClosingFilter3D( LabelImageType::Pointer & img, unsigned int ballsize )
{
  LabelImageType::Pointer imgDilate = BinaryDilateFilter3D( img, ballsize );

  return BinaryErodeFilter3D( imgDilate, ballsize );
}

// Function to convert a point from std::vector to itk::Point
// this also performs the RAS -> LPS conversion necessary
// from slicer -> ITK
itk::Point<float, 3> convertStdVectorToITKPoint(const std::vector<float> & vec)
{
	itk::Point<float, 3> p;

	// convert RAS to LPS
	p[0] = -vec[0];
	p[1] = -vec[1];
	p[2] = vec[2];
	return p;
}
