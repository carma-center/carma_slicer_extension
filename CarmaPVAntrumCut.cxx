/*
 *  CarmaPVAntrumCut.cxx
 */

#include "CarmaPVAntrumCutCLP.h"
#include "CarmaPVAntrumCut.h"

#include <vector>
#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>

#include <itkImageRegionConstIterator.h>

namespace
{
}

int main(int argc, char *argv[])
{
	PARSE_ARGS;
	
	// set up ITK types...
	typedef float InputPixelType;
	typedef unsigned int OutputPixelType;
	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	typedef itk::Image< OutputPixelType,  3 >   OutputImageType;
	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
	
	// read in the images...
	std::vector<InputImageType::Pointer> images;
	
	ReaderType::Pointer reader1 = ReaderType::New();
	reader1->SetFileName( endoLayer );
	InputImageType::Pointer image1 = reader1->GetOutput();
	reader1->Update();
	images.push_back(image1);
	
	ReaderType::Pointer reader2 = ReaderType::New();
	reader2->SetFileName( wallLayer );
	InputImageType::Pointer image2 = reader2->GetOutput();
	reader2->Update();
	images.push_back(image2);
	
	// Define the type of filter that we use.
	typedef itk::CarmaPVAntrumCut< InputImageType, InputImageType > filter_type;
	
	// Create a new ITK filter instantiation.	
	filter_type::Pointer filter = filter_type::New();
	
	filter->SetInput( 0, images[0] );
	filter->SetInput( 1, images[1] );
	
	filter->Update();
	
	// convert to mask
	InputImageType::Pointer output = filter->GetOutput();
	
	// create a new image to store the mask
	OutputImageType::Pointer mask_image = OutputImageType::New(); 
	mask_image->SetRegions( output->GetLargestPossibleRegion() ); 
	mask_image->Allocate();
	mask_image->SetOrigin( output->GetOrigin() );
	mask_image->SetSpacing( output->GetSpacing() );
	mask_image->FillBuffer ( itk::NumericTraits<OutputImageType::PixelType>::Zero ); 
	
	// iterate through data image, creating mask image
	ConstIteratorType it(output, output->GetLargestPossibleRegion());
	it.GoToBegin();
	for( ; !it.IsAtEnd(); ++it) {
		InputPixelType pixel = it.Get();
		if (pixel > 0) {
			mask_image->SetPixel(it.GetIndex(), 2);
		}
	}
	
	WriterType::Pointer writer;
	writer = WriterType::New();
	writer->SetFileName( endoNoVeins );
	writer->SetInput(mask_image);
	writer->UseCompressionOn();
	writer->Update();
	
	return EXIT_SUCCESS;
}

