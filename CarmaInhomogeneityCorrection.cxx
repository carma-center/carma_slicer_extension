//std
#include <string>


// auto-generated
#include "CarmaInhomogeneityCorrectionCLP.h"


#include "CarmaInhomogeneityCorrection.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// setup ITK types...
typedef float PixelType;
typedef unsigned int MaskPixelType;
typedef itk::Image< PixelType,  3 >   ImageType;
typedef itk::Image< MaskPixelType,  3 >   MaskImageType;
typedef itk::ImageFileReader< ImageType  >  ReaderType;
typedef itk::ImageFileWriter< ImageType  >  WriterType;

typedef InhomogeneityCorrectionFilter< ImageType, ImageType > CorrectionFilterType;


int main(int argc, char** argv)
{
  PARSE_ARGS; 

  // read data
  ReaderType::Pointer data_reader = ReaderType::New();
  data_reader->SetFileName(dataFileName);
  ImageType::Pointer data_image = data_reader->GetOutput();
  data_reader->Update();

  // read endo
  ReaderType::Pointer endo_reader = ReaderType::New();
  endo_reader->SetFileName(endoFileName);
  ImageType::Pointer endo_image = endo_reader->GetOutput();
  endo_reader->Update();

  // read roi
  ReaderType::Pointer roi_reader = ReaderType::New();
  roi_reader->SetFileName(roiFileName);
  ImageType::Pointer roi_image = roi_reader->GetOutput();
  roi_reader->Update();

  CorrectionFilterType::Pointer correction_filter = CorrectionFilterType::New();
  correction_filter->SetInput(0, data_image);
  correction_filter->SetInput(1, endo_image);
  correction_filter->SetInput(2, roi_image);

  //WriterType::Pointer writer;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFileName);
  writer->SetInput(correction_filter->GetOutput());
  writer->UseCompressionOn();
  writer->Update();

  return 0;
}
