/**
 * CarmaAutomaticLeftAtrialScar - automatic scar algorithm.
 * For details see paper: Perry, D. et al. "Automatic classification of scar tissue...".
 *                        Proceedings of SPIE Medical Imaging: Computer Aided Diagnosis. Feb 2012.
 */

//std
#include <string>

// auto-generated
#include "CarmaAutomaticLeftAtrialScarCLP.h"

// itk
#include <itkImage.h>
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <itkScalarImageKmeansImageFilter.h>

struct Cluster
{
  int count;
  double mean;
  Cluster(int count_=0, double mean_=0.f)
  :count(count_),
  mean(mean_)
  {}
};

int main(int argc, char** argv)
{

  PARSE_ARGS; 

  const int K = 4; 

  // setup ITK types...
  typedef float  PixelType;
  typedef itk::Image< PixelType,  3 >   ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;

  typedef unsigned char MaskPixelType;
  typedef itk::Image< MaskPixelType,  3 >   MaskImageType;
  typedef itk::ImageFileReader< MaskImageType  >  MaskReaderType;
 
  // read in the nrrds...
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( lgefn );
  ImageType::Pointer lge = reader->GetOutput();
  try
  {
    reader->Update();
  }
  catch(itk::ExceptionObject e)
  {
    std::cerr << "Error reading file " << lgefn << ": " << e << std::endl;
    exit(1);
  }
  ImageType::SizeType lgesize = lge->GetLargestPossibleRegion().GetSize();

  MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( lawallfn );
  MaskImageType::Pointer lawall = maskreader->GetOutput();
  try
  {
    maskreader->Update();
  }
  catch(itk::ExceptionObject e)
  {
    std::cerr << "Error reading file " << lawallfn << ": " << e << std::endl;
    exit(1);
  }
  MaskImageType::SizeType lawallsize = lawall->GetLargestPossibleRegion().GetSize();

  if( lawallsize != lgesize )
  {
    std::cerr << "Error: " << lgefn << " " << lgesize << " and " << lawallfn  << " " << lawallsize << " must have the same dimensions." << std::endl;
    exit(1);
  }

  maskreader = MaskReaderType::New();
  maskreader->SetFileName( laendofn );
  MaskImageType::Pointer laendo = maskreader->GetOutput();
  try
  {
    maskreader->Update();
  }
  catch(itk::ExceptionObject e)
  {
    std::cerr << "Error reading file " << laendofn << ": " << e << std::endl;
    exit(1);
  }
  MaskImageType::SizeType laendosize = laendo->GetLargestPossibleRegion().GetSize();

  if( lawallsize != laendosize )
  {
    std::cerr << "Error: " << laendofn << " " << laendosize << " and " << lawallfn  << " " << lawallsize << " must have the same dimensions." << std::endl;
    exit(1);
  }

  typedef itk::ImageRegionConstIterator< MaskImageType > ConstMaskRegionIterator;
  ConstMaskRegionIterator wallit( lawall, lawall->GetLargestPossibleRegion() );
  ConstMaskRegionIterator endoit( laendo, laendo->GetLargestPossibleRegion() );

  // calculate mean and standard deviation for statistical normalization:
  double mean = 0.f;
  double stdev = 0.f;
  int count = 0;

  for( wallit.GoToBegin(),endoit.GoToBegin(); !wallit.IsAtEnd(); ++wallit,++endoit )
  {
    if( wallit.Get() > 0 || endoit.Get() > 0) 
    {
      float p = lge->GetPixel(wallit.GetIndex());
      mean += p;
      ++count;
    }
  }
  mean /= count;

  for( wallit.GoToBegin(),endoit.GoToBegin(); !wallit.IsAtEnd(); ++wallit,++endoit )
  {
    if( wallit.Get() > 0 || endoit.Get() > 0) 
    {
      float p = lge->GetPixel(wallit.GetIndex());
      stdev += (p-mean)*(p-mean);
    }
  }
  stdev /= (count-1);
  stdev = std::sqrt( stdev );

  // we need to pull out the data in the wall+endo mask so that k-means only works on it
  // and not the 0 data outside the mask.. for speed and to match previous work.
  typedef itk::Image<float,1> ListType;
  typedef itk::ImageRegionIterator< ListType > ListIterator;
  ListType::Pointer epidata = ListType::New(); 
  ListType::Pointer wallmasklist = ListType::New(); 
  ListType::SizeType size;
  size.Fill(count);
  ListType::RegionType listregion(size);

  epidata->SetRegions( listregion ); 
  epidata->Allocate();
  ListIterator epidatait( epidata, epidata->GetLargestPossibleRegion() );

  wallmasklist->SetRegions( listregion ); 
  wallmasklist->Allocate();

  float min = 10000.f;
  float max = 0.f;
  for( wallit.GoToBegin(),endoit.GoToBegin(),epidatait.GoToBegin(); !wallit.IsAtEnd(); ++wallit,++endoit )
  {
    if( wallit.Get() > 0 || endoit.Get() > 0) 
    {
      float p = lge->GetPixel(wallit.GetIndex());
      p = (p-mean)/stdev; // statistically normalize the mri values:
      epidata->SetPixel(epidatait.GetIndex(), p );
      wallmasklist->SetPixel(epidatait.GetIndex(), wallit.Get());

      if( p > max ) max = p;
      if( p < min ) min = p;
      ++epidatait;
    }
  }

  // apply the itk k-means filter
  typedef itk::ScalarImageKmeansImageFilter< ListType, ListType > KmeansFilter;
  KmeansFilter::Pointer kmeansfilter = KmeansFilter::New();
  kmeansfilter->SetInput( epidata );
  // TODO: initialize randomly and run multiple times, then take result with best compactness like opencv.
  // see http://opencv.willowgarage.com/documentation/cpp/clustering_and_search_in_multi-dimensional_spaces.html
  // For now, initializing with evenly spaced values...
  float range = max - min;
  float subrange = range / 4; // sub1:[min,min+subrange), sub2:[min+subrange,min+2*subrange), ...
  for(int k=0; k<K; ++k)
  {
    kmeansfilter->AddClassWithInitialMean(min+k*subrange+(subrange/2)); // middle of each subrange, as defined above.
  }
  try
  {
    kmeansfilter->Update();
  }
  catch(itk::ExceptionObject e)
  {
    std::cerr << "Error running kmeans : " << e << std::endl;
    exit(1);
  }
  ListType::Pointer label = kmeansfilter->GetOutput();
 
  // analyze the labels, figure out which has the highest mean intensity:
  typedef std::map<int, Cluster> ClusterMap;
  ClusterMap clusterMap;
  int wallCount = 0;
  for( epidatait.GoToBegin(); !epidatait.IsAtEnd(); ++epidatait )
  {
    if( wallmasklist->GetPixel(epidatait.GetIndex()) > 0 )
    {
      ++wallCount;
      int l = static_cast<int>( label->GetPixel(epidatait.GetIndex()) );
      float p = epidata->GetPixel(epidatait.GetIndex());
      ClusterMap::iterator e = clusterMap.find(l);
      if( e == clusterMap.end() )
      {
        clusterMap[ l ].mean = p;
        clusterMap[ l ].count = 1;
      }
      else
      {
        clusterMap[ l ].mean += p;
        clusterMap[ l ].count += 1;
      }
    }
  }

  double maxVal = 0.f;
  int maxL = -1;
  for( ClusterMap::iterator e = clusterMap.begin(); e != clusterMap.end(); ++e )
  {
    (e->second).mean /= (e->second).count;
    if( (e->second).mean > maxVal )
    {
      maxVal =  (e->second).mean;
      maxL = (e->first);
    }
  }

  // send some results to std::cout:
  //std::cerr << "scar percentage: ";
  //std::cout << clusterMap[maxL].count / static_cast<float>(wallCount);
  //std::cerr << std::endl;

  // store the results in a nrrd: 
  MaskImageType::Pointer scarimage = MaskImageType::New(); 
  scarimage->SetRegions( lawall->GetLargestPossibleRegion() ); 
  scarimage->Allocate();
  scarimage->SetOrigin( lawall->GetOrigin() );
  scarimage->SetSpacing( lawall->GetSpacing() );

  for( wallit.GoToBegin(),endoit.GoToBegin(),epidatait.GoToBegin(); !wallit.IsAtEnd(); ++wallit,++endoit )
  {
    if( wallit.Get() > 0 && label->GetPixel(epidatait.GetIndex()) == maxL )
    {
      scarimage->SetPixel( wallit.GetIndex(), 1 );
    }
    else
    {
      scarimage->SetPixel( wallit.GetIndex(), 0 );
    }
    if( wallit.Get() > 0 || endoit.Get() > 0 ) ++epidatait;
  }

  typedef itk::ImageFileWriter< MaskImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( scarimage );
  writer->SetFileName( outputfn );
  writer->UseCompressionOn();

  try
  {
    writer->Update();
  }
  catch(itk::ExceptionObject e)
  {
    std::cerr << "Error writing file " << outputfn << ": " << e << std::endl;
  }
  
  return 0;
}
