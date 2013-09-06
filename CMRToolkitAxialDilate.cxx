//std
#include <string>


// auto-generated
#include "CMRToolkitAxialDilateCLP.h"
#include "CMRToolkitAxialDilate.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSliceBySliceImageFilter.h>
#include "itkPluginUtilities.h"
#include <itkBinaryDilateImageFilter.h>
#include "itkBinaryBallStructuringElement.h"


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
  typedef itk::Image< PixelType,  2 >   TwoDImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;

  // read mask
  typename ReaderType::Pointer mask_reader = ReaderType::New();
  mask_reader->SetFileName(targetFileName);
  typename ImageType::Pointer mask_image = mask_reader->GetOutput();
  mask_reader->Update();

  

  typedef itk::BinaryBallStructuringElement<
    PixelType,2>                  StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(4);
  structuringElement.CreateStructuringElement();


  typedef itk::BinaryDilateImageFilter< TwoDImageType, TwoDImageType, StructuringElementType> BinaryDilateImageFilterType;

  typedef itk::SliceBySliceImageFilter< ImageType, ImageType > SliceBySliceFilterType;

  typename BinaryDilateImageFilterType::Pointer dilate_filter = BinaryDilateImageFilterType::New();

  dilate_filter->SetKernel( structuringElement );
  dilate_filter->SetDilateValue( 1 );

  typename SliceBySliceFilterType::Pointer filter = SliceBySliceFilterType::New();

  filter->SetDimension( 2 );
  filter->SetFilter( dilate_filter );
  filter->SetInput( mask_image );

  /* 

  // 3D dilation

  typedef itk::BinaryBallStructuringElement<
    ImageType::PixelType,3>                  StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(2);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType> BinaryDilateImageFilterType;

  BinaryDilateImageFilterType::Pointer dilate_filter = BinaryDilateImageFilterType::New();
  dilate_filter->SetDilateValue( 1 );
  dilate_filter->SetKernel( structuringElement );
  dilate_filter->SetInput( mask_image );
  */

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFileName);
  //writer->SetInput(dilate_filter->GetOutput());
  writer->SetInput(filter->GetOutput());
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
    itk::GetImageType(targetFileName, pixelType, componentType);

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
