//std
#include <string>


// auto-generated
#include "CarmaInhomogeneityCorrectionCLP.h"
#include "CarmaInhomogeneityCorrection.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkPluginUtilities.h"



// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  typedef T PixelType;
  typedef itk::Image< PixelType,  3 >   ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;

  typedef InhomogeneityCorrectionFilter< ImageType, ImageType > CorrectionFilterType;


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

} // end of anonymous namespace


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
  {
    itk::GetImageType(dataFileName, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
    {
    case itk::ImageIOBase::UCHAR:
      return DoIt( argc, argv, static_cast<unsigned char>(0) );
      break;
    case itk::ImageIOBase::CHAR:
      return DoIt( argc, argv, static_cast<char>(0) );
      break;
    case itk::ImageIOBase::USHORT:
      return DoIt( argc, argv, static_cast<unsigned short>(0) );
      break;
    case itk::ImageIOBase::SHORT:
      return DoIt( argc, argv, static_cast<short>(0) );
      break;
    case itk::ImageIOBase::UINT:
      return DoIt( argc, argv, static_cast<unsigned int>(0) );
      break;
    case itk::ImageIOBase::INT:
      return DoIt( argc, argv, static_cast<int>(0) );
      break;
    case itk::ImageIOBase::ULONG:
      return DoIt( argc, argv, static_cast<unsigned long>(0) );
      break;
    case itk::ImageIOBase::LONG:
      return DoIt( argc, argv, static_cast<long>(0) );
      break;
    case itk::ImageIOBase::FLOAT:
      return DoIt( argc, argv, static_cast<float>(0) );
      break;
    case itk::ImageIOBase::DOUBLE:
      return DoIt( argc, argv, static_cast<double>(0) );
      break;
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
      std::cout << "unknown component type" << std::endl;
      break;
    }
  }

  catch( itk::ExceptionObject & excep )
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
