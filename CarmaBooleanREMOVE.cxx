/*
 * Boolean REMOVE filter to subtract two segmentations. 
 */

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkConstrainedValueDifferenceImageFilter.h"

#include "CarmaBooleanREMOVECLP.h"

int main( int argc, char * argv[] )
{
	PARSE_ARGS;

  typedef float InputPixelType;
  typedef float OutputPixelType;

  typedef itk::Image<InputPixelType, 3> InputImageType;
  typedef itk::Image<OutputPixelType, 3> OutputImageType;

  typedef itk::ImageFileReader<InputImageType> ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;

  typedef itk::BSplineInterpolateImageFunction<InputImageType> Interpolator;
  typedef itk::ResampleImageFilter<InputImageType, OutputImageType> ResampleType;
  typedef itk::ConstrainedValueDifferenceImageFilter<InputImageType, OutputImageType, OutputImageType> DifferenceFilterType;

  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
	
  reader1->SetFileName( inputVolume1 );
  reader2->SetFileName( inputVolume2 );

  reader1->Update();
  reader2->Update();

  Interpolator::Pointer interp = Interpolator::New();
  interp->SetInputImage( reader2->GetOutput() );
  interp->SetSplineOrder( 1 );

  ResampleType::Pointer resample = ResampleType::New();
  resample->SetInput( reader2->GetOutput() );
  resample->SetOutputParametersFromImage( reader1->GetOutput() );
  resample->SetInterpolator( interp );
  resample->SetDefaultPixelValue( 0 );
	resample->Update();
 
  DifferenceFilterType::Pointer differenceFilter = DifferenceFilterType::New();
  differenceFilter->SetInput1( reader1->GetOutput() );
  differenceFilter->SetInput2( resample->GetOutput() );
	differenceFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume );
  writer->SetInput( differenceFilter->GetOutput() );
  writer->Update();
		
	return EXIT_SUCCESS;
}
