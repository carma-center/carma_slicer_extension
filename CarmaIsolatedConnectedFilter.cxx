/*
 *  CarmaIsolatedConnectedFilter.cxx
 */

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkIsolatedConnectedImageFilter.h"
#include "CarmaIsolatedConnectedFilterCLP.h"
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

int main( int argc, char *argv[] )
{
  PARSE_ARGS;
	
	typedef unsigned char                            OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastingFilterType;
  
  CastingFilterType::Pointer caster = CastingFilterType::New();
  typedef  itk::ImageFileReader< InternalImageType > ReaderType;
  typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( inputImage );
  writer->SetFileName( outputImage );
	InternalImageType::Pointer inputImagePtr = InternalImageType::New();
	inputImagePtr = reader->GetOutput();
  typedef itk::CurvatureFlowImageFilter< InternalImageType, InternalImageType > CurvatureFlowImageFilterType;
  CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();

  typedef itk::IsolatedConnectedImageFilter< InternalImageType, InternalImageType > ConnectedFilterType;

  ConnectedFilterType::Pointer isolatedConnected = ConnectedFilterType::New();
  reader->Update();
  smoothing->SetInput( reader->GetOutput() );
  isolatedConnected->SetInput( smoothing->GetOutput() );
  caster->SetInput( isolatedConnected->GetOutput() );
  writer->SetInput( caster->GetOutput() );
  smoothing->SetNumberOfIterations( 5 );
  smoothing->SetTimeStep( 0.125 );

	PointList seeds1_list;
	PointList seeds2_list;
	
	typedef  std::vector<InternalImageType::IndexType> IndexList;
	
  isolatedConnected->FindUpperThresholdOn();
	isolatedConnected->SetLower( lowerThreshold );
	
	IndexList indexList1;
	IndexList indexList2;

  if( seeds1.size() > 0 && seeds2.size() > 0 )
	{
	  seeds1_list.resize( seeds1.size() );
		seeds2_list.resize( seeds2.size() );
		
		// Convert both point lists to ITK points and convert RAS -> LPS
		std::transform( seeds1.begin(), seeds1.end(),
                   seeds1_list.begin(), convertStdVectorToITKPoint );
									
    std::transform( seeds2.begin(), seeds2.end(),
                   seeds2_list.begin(), convertStdVectorToITKPoint );
								
	  indexList1.resize( seeds1.size() );
		indexList2.resize( seeds2.size() );
		
		for( unsigned int i = 0; i<seeds1.size(); i++ )
		{
			inputImagePtr->TransformPhysicalPointToIndex( seeds1_list[i], indexList1[i] );
			isolatedConnected->SetSeed1( indexList1[i] );
		}
		
		for( unsigned int i = 0; i<seeds2.size(); i++ )
		{
			inputImagePtr->TransformPhysicalPointToIndex( seeds2_list[i], indexList2[i] );
			isolatedConnected->SetSeed2( indexList2[i] );
		}							
	}
  else
	{
	  std::cerr << "Must supply at least one seed for each region" << std::endl;
    return EXIT_FAILURE;
	}

  isolatedConnected->SetReplaceValue( 255 );
	
	if( inputImagePtr->GetBufferedRegion().IsInside( indexList1[0] ) )
	{
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
	}
	else
	{
		std::cerr << "Index not within image boundaries" << std::endl;
	}
	
	return EXIT_SUCCESS;
}