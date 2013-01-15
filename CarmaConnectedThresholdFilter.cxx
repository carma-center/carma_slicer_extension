/*
 *  CarmaConnectedThresholdFilter.cxx
 */

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "CarmaConnectedThresholdFilterCLP.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

namespace
{
  typedef   float           InternalPixelType;
  const     unsigned int    Dimension = 3;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
	typedef InternalImageType::PointType ImagePointType;
	typedef  std::vector<ImagePointType>        PointList;

  // Method to convert RAS point to LPS point
  static itk::Point<float, 3> convertStdVectorToITKPoint(const std::vector<float> & vec)
  {
    itk::Point<float, 3> p;

    // convert RAS to LPS
    p[0] = -vec[0];
    p[1] = -vec[1];
    p[2] = vec[2];
    return p;
  }
}

int main( int argc, char *argv[])
{
  PARSE_ARGS;

  const  unsigned int  Dimension = 3;
	typedef  float  InternalPixelType;
	typedef	 unsigned char  OutputPixelType;
  typedef	 itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef	 itk::Image< OutputPixelType, Dimension >  OutputImageType;
  typedef	 itk::CastImageFilter< InternalImageType, OutputImageType >  CastingFilterType;
  typedef	 itk::CurvatureFlowImageFilter< InternalImageType, InternalImageType >  CurvatureFlowImageFilterType;
	typedef	 itk::ConnectedThresholdImageFilter< InternalImageType, InternalImageType >  ConnectedFilterType;
  typedef  itk::ImageFileReader< InternalImageType >  ReaderType;
  typedef  itk::ImageFileWriter<  OutputImageType  >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

	reader->SetFileName( inputImage );
	InternalImageType::Pointer inputImagePtr = InternalImageType::New();
	inputImagePtr = reader->GetOutput();
	reader->Update();
	writer->SetFileName( outputImage );
	
  CastingFilterType::Pointer caster = CastingFilterType::New();
  ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
  CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();
  smoothing->SetInput( inputImagePtr );
  connectedThreshold->SetInput( smoothing->GetOutput() );
  caster->SetInput( connectedThreshold->GetOutput() );
  writer->SetInput( caster->GetOutput() );

  smoothing->SetNumberOfIterations( 5 );
  smoothing->SetTimeStep( 0.125 );

  connectedThreshold->SetLower(  lowerThreshold  );
	connectedThreshold->SetUpper(  upperThreshold  );	
  connectedThreshold->SetReplaceValue( 255 );
	
	PointList seeds_list;
	
	typedef  std::vector<InternalImageType::IndexType> IndexList;
	IndexList indexList;
	
	if( seeds.size() > 0 )
	{
	  seeds_list.resize( seeds.size() );
		
		// Convert both point lists to ITK points and convert RAS -> LPS
		std::transform( seeds.begin(), seeds.end(), seeds_list.begin(), convertStdVectorToITKPoint );
									 
		indexList.resize( seeds.size() );
	  
		for( unsigned int i = 0; i<seeds.size(); i++)
		{
			bool isInside = inputImagePtr->TransformPhysicalPointToIndex( seeds_list[i], indexList[i] );
			if(isInside)
			{
				connectedThreshold->SetSeed( indexList[i] );
			}
			else
			{
			  std::cout << "Seed is outside image region." << std::endl;
				return EXIT_FAILURE;
			}	
		}
	}
	else
	{
		std::cerr << "Must supply at least one seed for segmentation." << std::endl;
		return EXIT_FAILURE;
  }
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
  }

  return EXIT_SUCCESS;
}
