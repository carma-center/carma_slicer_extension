
#ifndef _itkLeftAtriumIsolationFilter_H
#define _itkLeftAtriumIsolationFilter_H

// std includes
#include <algorithm>

// itk includes
#include <itkConfigure.h>
#include <itkNumericTraits.h>
#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include <itkImageRegionIterator.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkZeroFluxNeumannBoundaryCondition.h>
#include <itkOffset.h>
#include <itkProgressReporter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkSliceBySliceImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkCommand.h>
 
namespace itk
{
  
  using itk::ImageToImageFilter;
  using itk::SmartPointer;

  class ProgressObserver : public itk::Command
  {
  public:
    ProgressObserver( itk::ProgressReporter & rep )
    :reporter(rep)
    {}
    
    virtual ~ProgressObserver() {}

    virtual void Execute( itk::Object *caller, const itk::EventObject & event )
    {
      reporter.CompletedPixel();
    }

    virtual void Execute(const itk::Object *caller, const itk::EventObject & event )
    {
      reporter.CompletedPixel();
    }

  private:
    itk::ProgressReporter & reporter;
  };

  template <class TInputImage, class TOutputImage=TInputImage>
  class itkLeftAtriumIsolationFilter :
    public ImageToImageFilter < TInputImage, TOutputImage >
  {
    public:
      /** Extract dimension from input and output image. */
      itkStaticConstMacro(InputImageDimension, unsigned int,
                          TInputImage::ImageDimension);
      itkStaticConstMacro(OutputImageDimension, unsigned int,
                          TOutputImage::ImageDimension);

      /** Convenient typedefs for simplifying declarations. */
      typedef TInputImage  InputImageType;
      typedef TOutputImage OutputImageType;

      /** Standard class typedefs. */
      typedef itkLeftAtriumIsolationFilter Self;
      typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
      typedef SmartPointer<Self> Pointer;
      typedef SmartPointer<const Self> ConstPointer;

      /** Method for creation through the object factory. */
      itkNewMacro(Self);

      /** Run-time type information (and related methods). */
      itkTypeMacro(itkLeftAtriumIsolationFilter, ImageToImageFilter);
  
      /** Image typedef support. */
      typedef typename InputImageType::PixelType  InputPixelType;
      typedef typename OutputImageType::PixelType OutputPixelType;
    
      typedef typename InputImageType::RegionType  InputImageRegionType;
      typedef typename OutputImageType::RegionType OutputImageRegionType;

      typedef typename InputImageType::SizeType InputSizeType;

      // get/set-ers:
      itkSetMacro(NoValue, InputPixelType);
      itkGetConstReferenceMacro(NoValue, InputPixelType);
      itkSetMacro(InterestPointValue, OutputPixelType);
      itkGetConstReferenceMacro(InterestPointValue, OutputPixelType);
      itkSetMacro(SliceIndex, unsigned int);
      itkGetConstReferenceMacro(SliceIndex, unsigned int);
      itkSetMacro(FilterRadius, InputSizeType);
      //itkGetConstReferenceMacro(FilterRadius, InputSizeType);
      InputSizeType GetFilterRadius() const
      {
        InputSizeType r = m_FilterRadius;
        r[m_SliceIndex] = 0;
        return r;
      }

     /** itkLeftAtriumIsolationFilter needs a larger input requested region than
       * the output requested region.  As such, itkLeftAtriumIsolationFilter needs
       * to provide an implementation for GenerateInputRequestedRegion()
       * in order to inform the pipeline execution model.
       *
       * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
      //virtual void GenerateInputRequestedRegion() throw(itk::InvalidRequestedRegionError);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck,
                  (itk::Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  itkConceptMacro(InputConvertibleToOutputCheck,
                  (itk::Concept::Convertible<InputPixelType, OutputPixelType>));
  itkConceptMacro(InputLessThanComparableCheck,
                  (itk::Concept::LessThanComparable<InputPixelType>));
  /** End concept checking */
#endif

 protected:
    itkLeftAtriumIsolationFilter();
    virtual ~itkLeftAtriumIsolationFilter() {}

    void PrintSelf(std::ostream &os, itk::Indent indent) const;

    void GenerateData();

  private:
    itkLeftAtriumIsolationFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    InputSizeType m_FilterRadius;
    InputPixelType m_NoValue;
    OutputPixelType m_InterestPointValue;
    unsigned int m_SliceIndex; // which dimension will iterate over slices
 }; // end class itkLeftAtriumIsolationFilter

  // DEFINITIONS:

  template <class TInputImage, class TOutputImage>
  itkLeftAtriumIsolationFilter<TInputImage, TOutputImage>
  ::itkLeftAtriumIsolationFilter()
  :m_FilterRadius(),
  m_NoValue(itk::NumericTraits<InputPixelType>::Zero),
  m_InterestPointValue(itk::NumericTraits<OutputPixelType>::Zero + 1),
  m_SliceIndex(2) // default slice dimension is "z"
  {
    m_FilterRadius.Fill(1);
  }

  // callback for penalizing a score if the line connecting two points crosses the wall mask
  template< class ImageType >
  void penalizeMaskCrossing( typename ImageType::ConstPointer & mask, const typename ImageType::IndexType & pt, float & score , float penalty)
  {
    if( mask->GetPixel(pt) > 0 )
    {
      score += penalty;
    }
  }


  // callback for setting the pixel value in bresenham(..) below.
  template< class OutputImageType >
  void setPixel( typename OutputImageType::Pointer & im, const typename OutputImageType::IndexType & in, const typename OutputImageType::PixelType & p )
  {
    im->SetPixel( in, p );
  }

  /** Bresenham's algorithm for drawing a line in a slice
   * @param pair the pair of points making up the endpoints to the line.
   * @param outputCb function called for each point of the line.
   * NOTE: assumes SliceIndex = 2, TODO: make it work for any slice index...
   */
  template< class IndexType, class SizeType, class ImageType >
  void bresenham( const IndexType & start, const IndexType & end, const SizeType & imageSize, typename ImageType::ConstPointer & mask, float & score, float penalty )
  {
      int p0[3], p1[3];
			
      p0[0] = static_cast<int>( start[0] );
      p0[1] = static_cast<int>( start[1] );
      p0[2] = static_cast<int>( start[2] );
      
      p1[0] = static_cast<int>( end[0] );
      p1[1] = static_cast<int>( end[1] );
      p1[2] = static_cast<int>( end[2] );

      int dx = std::abs( p1[0] - p0[0] );
      int dy = std::abs( p1[1] - p0[1] );
      int sx,sy;
      if ( p0[0] < p1[0] )
      {
        sx = 1;
      }
      else
      {
        sx = -1;
      }
      if ( p0[1] < p1[1] )
      {
        sy = 1;
      }
      else
      {
        sy = -1;
      }

      int err = dx - dy;
      
      IndexType loc;
      loc[2] = start[2]; // should be on the same slice
      while ( true ) // see break statement below
      {
        loc[0] = static_cast<typename IndexType::IndexValueType>(p0[0]);
        loc[1] = static_cast<typename IndexType::IndexValueType>(p0[1]);
        if ( loc[0] >= imageSize[0] || loc[1] >= imageSize[1] || p0[0] < 0.f || p0[1] < 0.f )
        {
          // reached an invalid range... stop.
          break;
        }
				penalizeMaskCrossing<ImageType>(mask, loc, score, penalty);
				
        if ( ((sx == 1 && p0[0] >= p1[0]) || (sx == -1 && p0[0] <= p1[0])) && 
             ((sy == 1 && p0[1] >= p1[1]) || (sy == -1 && p0[1] <= p1[1])) ) 
        {
          // completed the line drawing... stop.
          break;
        }
        int e2 = err * 2;
        if ( e2 > -dy )
        {
          err -= dy;
          p0[0] += sx;
        }
        if ( e2 < dx )
        {
          err += dx;
          p0[1] += sy;
        }
      }

  }

  /*template< class IndexType , class SizeType, class PixelType>//, class ImageType >
  void bresenham2(const IndexType & start, const IndexType & end, const SizeType & imageSize, float p)//, typename ImageType::ConstPointer & im, const typename ImageType::PixelType & p)
  {
      int p0[3], p1[3];
			
      p0[0] = static_cast<int>( start[0] );
      p0[1] = static_cast<int>( start[1] );
      p0[2] = static_cast<int>( start[2] );
      //float p = 20;
      p1[0] = static_cast<int>( end[0] );
      p1[1] = static_cast<int>( end[1] );
      p1[2] = static_cast<int>( end[2] );

      int dx = std::abs( p1[0] - p0[0] );
      int dy = std::abs( p1[1] - p0[1] );
      int sx,sy;
      if ( p0[0] < p1[0] )
      {
        sx = 1;
      }
      else
      {
        sx = -1;
      }
      if ( p0[1] < p1[1] )
      {
        sy = 1;
      }
      else
      {
        sy = -1;
      }

      int err = dx - dy;
      
      IndexType loc;
      loc[2] = start[2]; // should be on the same slice
      while ( true ) // see break statement below
      {
        loc[0] = static_cast<typename IndexType::IndexValueType>(p0[0]);
        loc[1] = static_cast<typename IndexType::IndexValueType>(p0[1]);
        if (loc[0] >= imageSize[0] || loc[1] >= imageSize[1] || p0[0] < 0.f || p0[1] < 0.f)
        {
          // reached an invalid range... stop.
          break;
        }
        //outputCb(loc);
				//setPixel<ImageType>(im, loc, p);
				
        if ( ((sx == 1 && p0[0] >= p1[0]) || (sx == -1 && p0[0] <= p1[0])) && 
             ((sy == 1 && p0[1] >= p1[1]) || (sy == -1 && p0[1] <= p1[1])) ) 
        {
          // completed the line drawing... stop.
          break;
        }
        int e2 = err * 2;
        if (e2 > -dy)
        {
          err -= dy;
          p0[0] += sx;
        }
        if (e2 < dx )
        {
          err += dx;
          p0[1] += sy;
        }
      }

  }*/


  /**
   * Calculate euclidean distance from an itk index type
   */
  template< class IndexType >
  float dist(const IndexType & p1, const IndexType & p2)
  {
    float d0 = static_cast<float>( p1[0] - p2[0] );
    float d1 = static_cast<float>( p1[1] - p2[1] );
    float d2 = static_cast<float>( p1[2] - p2[2] );
    return std::sqrt( d0*d0 + d1*d1 + d2*d2 );
  }

  /**
   * Add a candidate to the list, putting in the list with other similar (clustered by Euclidean distance) points.
   */
  template< class ElementType, class ListType >
  void addCandidate( const ElementType & e, std::vector<ListType> & list, unsigned int sliceIndex)
  {
    for( unsigned int i=0; i< list.size(); ++i )
    {
      if( e[sliceIndex] == list[i][0][sliceIndex] && dist( e, list[i][0] ) < 10.0f ) 
      {
        // probably a member of this cluster list
        // now verify:
        for( unsigned int j=0; j < list[i].size(); ++j  )
        {
          if( dist( e, list[i][j] ) < 3.0f )
          {
            list[i].push_back( e );
            return;
          }
        }
      }
    }
    
    // didn't find any matches, must be a new cluster...
    ListType newList;
    newList.push_back( e );
    list.push_back( newList );
  }

  /**
   * reduce a list of candidate lists, to a list of averages of those candidates.
   */
  template< class ListType >
  ListType consolidateCandidates( const std::vector<ListType> & candidateList )
  {
    typedef typename ListType::value_type ElementType;
    typedef typename ElementType::IndexValueType NumberType;

    ListType finalList;
    for( unsigned int i=0; i<candidateList.size(); ++i )
    {
      ElementType mean;
      float xSum = 0.f;
      float ySum = 0.f;
      float zSum = 0.f;
      unsigned int len = candidateList[i].size();
      for ( unsigned int j=0; j<len; ++j )
      {
        xSum += static_cast<float>( candidateList[i][j][0] );
        ySum += static_cast<float>( candidateList[i][j][1] );
        zSum += static_cast<float>( candidateList[i][j][2] );
      }
      mean[0] = static_cast<NumberType>( xSum / len );
      mean[1] = static_cast<NumberType>( ySum / len );
      mean[2] = static_cast<NumberType>( zSum / len );
      finalList.push_back(mean);
    }
    return finalList;
  }

  /**
   * Nearest Neighbor pairing
   * @param list - list of points to be paired off.
   * @return list of point pairs.
   */
  template< class ListType, class ImageType >
  std::vector<ListType> pairOff( const ListType & list, unsigned int sliceIndex, typename ImageType::ConstPointer & wall )
  {
    typedef typename ListType::value_type IndexType;
    // first sort into lists per slice:
    std::vector<ListType> slices;
    typename std::vector<ListType>::iterator sliceIt;
    typename ListType::const_iterator outer;
    typename ListType::const_iterator inner;
    for ( outer = list.begin(); outer != list.end(); ++outer )
    {
      bool newSlice = true;
      for ( sliceIt = slices.begin(); sliceIt != slices.end() ; ++sliceIt )
      {
        if ( (*outer)[sliceIndex] == (*sliceIt)[0][sliceIndex] )
        {
          newSlice = false;
          break;
        }
      }
      if ( newSlice )
      {
        ListType slice;
        slice.push_back(*outer);
        inner = outer;
        for ( ++inner; inner != list.end(); ++inner )
        {
          if ( (*outer)[sliceIndex] == (*inner)[sliceIndex] )
          {
            slice.push_back(*inner);
          }
        }
        slices.push_back(slice);
      }
    }
    
    std::vector<ListType> pairList;
    typename ImageType::SizeType imageSize = wall->GetLargestPossibleRegion().GetSize();
    // brute force solve minimum overall distances...
    // ultimately we want pairList = ( (p1,p2), (p3,p4), ... ), 
    // such that dist(p1,p2) + dist(p3,p4) + ... is the smallest possible
    // (for any pairing).
    for ( sliceIt = slices.begin(); sliceIt != slices.end(); ++sliceIt )
    {
      unsigned int numPts = (*sliceIt).size();

      // now figure out the min pairList for this slice...
      if ( numPts == 2 )
      {
        ListType pair;
        pair.push_back( (*sliceIt)[0] );
        pair.push_back( (*sliceIt)[1] );        
        pairList.push_back(pair);
        continue;
      }

      // here is the brute force search... (we don't expect to ever have more than 10 points)...
      std::vector<unsigned int> indices(numPts);
      std::vector<unsigned int> bestIndices(numPts);
      float bestScore = 10000.f;
      for ( unsigned int i=0; i<numPts; ++i ) 
      {
        indices[i] = i;
        bestIndices[i] = i;
      }
      while ( std::next_permutation( indices.begin() , indices.end() ) )
      {
        float score = 0.f;
        for ( unsigned int i=0; i+1<numPts; i+=2 )
        {
          IndexType p0 = (*sliceIt)[indices[i]];
          IndexType p1 = (*sliceIt)[indices[i+1]];
          if ( i+1 >= numPts ) break; // odd number of pts...
          score += dist( p0 , p1 );

          // Penalty if any of the pairs cross the mask..
          float penalty = 10.f;
          bresenham<IndexType, typename ImageType::SizeType, ImageType>( p0, p1, imageSize, wall, score, penalty );
        }
        if (score < bestScore)
        {
          bestScore = score;
          for ( unsigned int i=0; i<numPts; ++i ) bestIndices[i] = indices[i];
        }
      }

      // translate into pairs:
      for ( unsigned int i=0; i+1<numPts; i+=2 )
      {
        ListType pair;
        pair.push_back( (*sliceIt)[ bestIndices[i] ] );
        pair.push_back( (*sliceIt)[ bestIndices[i+1] ] );
        pairList.push_back( pair );
      }
    }

    return pairList;
  }

  /**
   * Calc index
   */
  inline unsigned int ind(unsigned int r, unsigned int c, unsigned int width)
  {
    return (r*width) + c;
  }

  /**
   * Check if the mask data forms a line in the region.
   * @param [in] imageData - neighborhood iterator
   * @return true if a line is found in the region, false otherwise.
   */
  template< class NItType >
  bool lineTest( NItType region )
  {
    typename NItType::RadiusType radius = region.GetRadius();
    unsigned int width = 1 + 2 * radius[0];
    unsigned int height = 1 + 2 * radius[1];

    // check for row lines:
    for( unsigned int r=0; r<height; ++r )
    {
      bool line = true;
      for( unsigned int c=0; c<width; ++c )
      {
        if ( region.GetPixel(ind(r,c,width)) == 0 )
        {
          line = false;
          break;
        }
      }
      if (line) return true;
    }
    // check for col lines:
    for( unsigned int c=0; c<width; ++c )
    {
      bool line = true;
      for( unsigned int r=0; r<height; ++r )
      {
        if ( region.GetPixel(ind(r,c,width)) == 0 )
        {
          line = false;
          break;
        }
      }
      if (line) return true;
    }
    // check for top-bottom diagonal lines:
    for( unsigned int row=0; row<(height-1); ++row)
    {
      bool line = true;
      for( unsigned int r=row, c=0; c<width && r<height; ++c,++r )
      {
        if ( region.GetPixel(ind(r,c,width)) == 0 )
        {
          line = false;
          break;
        }
      }
      if (line) return true;
    }
    // check for bottom-top diagonal lines:
    for( unsigned int row=1; row<height; ++row)
    {
      bool line = true;
      for( int r=row, c=0; c<width && r>=0; ++c,--r )
      {
        if ( region.GetPixel(ind(r,c,width)) == 0 )
        {
          line = false;
          break;
        }
      }
      if (line) return true;
    }

    // no lines detected
    return false;
  }

  /**
   * Checks that both endo and wall neighborhoods have lines.
   */
  template< class WallNItType , class EndoNItType>
  bool lineTest( WallNItType wallRegion , EndoNItType endoRegion )
  {
    return lineTest(wallRegion) && lineTest(endoRegion);
  }

  /**
   * Checks if there is a 3-way boundary between wall, endo, and outer.
   */
  template< class WallNItType , class EndoNItType, class MaskImageType >
  bool triBoundaryTest( WallNItType wallRegion, MaskImageType wallImage, EndoNItType endoRegion, MaskImageType endoImage )
  {
    // center pixel must be a wall pixel
    if( wallRegion.GetCenterPixel() == 0.f )
    {
      return false;
    }

    typedef typename WallNItType::RadiusType RadiusType;
    typedef typename WallNItType::IndexType IndexType;
    RadiusType radius = wallRegion.GetRadius();
    IndexType center = wallRegion.GetIndex( wallRegion.GetCenterNeighborhoodIndex() );

    size_t count = (1 + 2 * radius[0]) * (1 + 2 * radius[1]);

    for(size_t i=0; i<count; ++i)
    {
      IndexType loc = endoRegion.GetIndex(i);
      if( loc == center )
      {
        continue; // skip center
      }
      if( endoRegion.GetPixel(i) > 0.f )
      {
        // check if this endo pixel borders any outer pixels (aka non-endo and non-wall)
        IndexType check = loc;
        //check north
        ++(check[0]);
        if( wallImage->GetPixel( check ) == 0.f && endoImage->GetPixel( check ) == 0.f ) return true;
        //check south
        check[0]-=2;
        if( wallImage->GetPixel( check ) == 0.f && endoImage->GetPixel( check ) == 0.f ) return true;
        //check east
        ++(check[0]);
        ++(check[1]);
        if( wallImage->GetPixel( check ) == 0.f && endoImage->GetPixel( check ) == 0.f ) return true;
        //check west
        check[1]-=2;
        if( wallImage->GetPixel( check ) == 0.f && endoImage->GetPixel( check ) == 0.f ) return true;
      }
    }
    return false;
  }

  template< class TInputImage, class TOutputImage>
  void
  itkLeftAtriumIsolationFilter< TInputImage, TOutputImage>
  ::GenerateData()
  {
    // Allocate output
    typename OutputImageType::Pointer output = this->GetOutput();
    output->SetBufferedRegion( output->GetLargestPossibleRegion() );   
    output->Allocate();
    output->FillBuffer ( itk::NumericTraits<OutputPixelType>::Zero ); 

    // Get Input pointers
    typename  InputImageType::ConstPointer inputEndo = this->GetInput(0);
    typename  InputImageType::ConstPointer inputWall = this->GetInput(1);

    // Find the data-set boundary "faces"
    typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> FacesCalculatorType;
    FacesCalculatorType bC;
    // assume all 3 images are same size...
    typename FacesCalculatorType::FaceListType faceList = bC( output, inputEndo->GetLargestPossibleRegion(), GetFilterRadius() );

    // support progress methods/callbacks
    unsigned int progressSectionLength = inputEndo->GetLargestPossibleRegion().GetNumberOfPixels();
    itk::ProgressReporter progress(this, 0, progressSectionLength * 3);

    // All of our neighborhoods have an odd number of pixels, so there is
    // always a median index (if there where an even number of pixels
    // in the neighborhood we have to average the middle two values).

    itk::ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;
    typedef itk::ConstNeighborhoodIterator<InputImageType> NeighborhoodIteratorType;
    typedef typename NeighborhoodIteratorType::IndexType IndexType;
    typedef std::vector<IndexType> IndexListType;
    typedef std::vector<IndexListType> ListOfIndexListType;

    ListOfIndexListType interestPointCandidates;

    //**** Part 1 - scan for important points

    // Process each of the boundary faces.  These are N-d regions which border
    // the edge of the buffer.
    for ( typename FacesCalculatorType::FaceListType::iterator fit=faceList.begin(); 
          fit != faceList.end(); ++fit)
    {
      NeighborhoodIteratorType endoIt(GetFilterRadius(), inputEndo, *fit);
      NeighborhoodIteratorType wallIt(GetFilterRadius(), inputWall, *fit);
      endoIt.OverrideBoundaryCondition(&nbc);
      wallIt.OverrideBoundaryCondition(&nbc);

      endoIt.GoToBegin();
      wallIt.GoToBegin();
      while ( ! (endoIt.IsAtEnd() || wallIt.IsAtEnd()) )
      {
        // setting output to wall values:
        output->SetPixel( wallIt.GetIndex() , wallIt.GetCenterPixel() );
        // Check all the pixels in the neighborhood, note that we use
        // GetPixel on the NeighborhoodIterator to honor the boundary conditions
        bool endoPixel = false;
        bool wallPixel = false;
        bool emptyPixel = false;
        unsigned int neighborhoodSize = endoIt.Size();
        IndexType wallLocation;
        for ( unsigned int i=0; i<neighborhoodSize; ++i )
        {
          InputPixelType wallValue = wallIt.GetPixel(i);
          InputPixelType endoValue = endoIt.GetPixel(i);

          if ( wallValue == m_NoValue && endoValue == m_NoValue )
          {
            emptyPixel = true;
          }
          else if ( wallValue != m_NoValue && endoValue != m_NoValue )
          {
            wallLocation = wallIt.GetIndex(i);
            wallPixel = true; // when overlapping count as wall pixel.
          }
          else if ( wallValue != m_NoValue )
          {
            wallLocation = wallIt.GetIndex(i);
            wallPixel = true;
          }
          else
          {
            endoPixel = true;
          }
        }
        if ( endoPixel && wallPixel && emptyPixel )
        {
          // checks that there is a wall, endo, and empty pixel next to each other:
          // (also center pixel has to be a wall pixel)
          if ( triBoundaryTest( wallIt, inputWall, endoIt, inputEndo ) )
          {
            addCandidate( wallIt.GetIndex(), interestPointCandidates, GetSliceIndex() );
          }
        }
        ++endoIt;
        ++wallIt;
        progress.CompletedPixel();
      }
    }
   

    //***** Part 2 - apply cutoffs to endo layer
    // we accomplish this by hole filling the union of lines and the wall

    // consolidate interest points, and pair them off:
    IndexListType interestPoints = consolidateCandidates( interestPointCandidates );
    ListOfIndexListType interestPairs = pairOff<IndexListType,InputImageType>( interestPoints, GetSliceIndex(), inputWall );
		typename OutputImageType::SizeType imageSize = output->GetLargestPossibleRegion().GetSize();
		
    for ( unsigned int i=0; i<interestPairs.size(); ++i )
    {
			typename OutputImageType::PixelType interestPointVal = this->GetInterestPointValue();
      IndexListType pair = interestPairs[i];
			int p0[3], p1[3];
			
      p0[0] = static_cast<int>( pair[0][0] );
      p0[1] = static_cast<int>( pair[0][1] );
      p0[2] = static_cast<int>( pair[0][2] );
      
      p1[0] = static_cast<int>( pair[1][0] );
      p1[1] = static_cast<int>( pair[1][1] );
      p1[2] = static_cast<int>( pair[1][2] );

      int dx = std::abs( p1[0] - p0[0] );
      int dy = std::abs( p1[1] - p0[1] );
      int sx,sy;
      if ( p0[0] < p1[0] )
      {
        sx = 1;
      }
      else
      {
        sx = -1;
      }
      if ( p0[1] < p1[1] )
      {
        sy = 1;
      }
      else
      {
        sy = -1;
      }

      int err = dx - dy;
      
      IndexType loc;
      loc[2] = pair[0][2]; // should be on the same slice
      while ( true ) // see break statement below
      {
        loc[0] = static_cast<typename IndexType::IndexValueType>(p0[0]);
        loc[1] = static_cast<typename IndexType::IndexValueType>(p0[1]);
        if (loc[0] >= imageSize[0] || loc[1] >= imageSize[1] || p0[0] < 0.f || p0[1] < 0.f)
        {
          // reached an invalid range... stop.
          break;
        }
				setPixel<OutputImageType>(output, loc, interestPointVal);
				
        if ( ((sx == 1 && p0[0] >= p1[0]) || (sx == -1 && p0[0] <= p1[0])) && 
             ((sy == 1 && p0[1] >= p1[1]) || (sy == -1 && p0[1] <= p1[1])) ) 
        {
          // completed the line drawing... stop.
          break;
        }
        int e2 = err * 2;
        if ( e2 > -dy )
        {
          err -= dy;
          p0[0] += sx;
        }
        if ( e2 < dx )
        {
          err += dx;
          p0[1] += sy;
        }
      }
    }

    //DEBUG:
    //return;
    
    // Flood fill the background with the foreground value (2D, slice by slice)
    typedef typename itk::Image<typename OutputImageType::PixelType , OutputImageType::ImageDimension-1> SliceImageType;
    typedef typename itk::ConnectedThresholdImageFilter<SliceImageType, SliceImageType> FloodFillType;
    typename FloodFillType::Pointer floodFilter = FloodFillType::New();
    typename FloodFillType::IndexType seed;
    seed.Fill(0); // corner pixel

    floodFilter->SetLower( GetNoValue() );
    floodFilter->SetUpper( GetNoValue() );
    floodFilter->SetSeed( seed );
    floodFilter->SetReplaceValue( GetInterestPointValue() );

    typedef itk::SliceBySliceImageFilter<OutputImageType, OutputImageType, FloodFillType, FloodFillType, SliceImageType, SliceImageType> SliceBySliceType;

    typename SliceBySliceType::Pointer sliceFilter = SliceBySliceType::New();

    sliceFilter->SetFilter( floodFilter );
    sliceFilter->SetInput( output );
    sliceFilter->AddObserver( itk::ProgressEvent(), new ProgressObserver(progress) );
    sliceFilter->Update();

    // Now set the result to be the inverse of the union of flood fill result and wall (and atrium closing lines)
    typename OutputImageType::Pointer result = sliceFilter->GetOutput();
    itk::ImageRegionIterator< OutputImageType > resultIt( result, result->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< OutputImageType > outputIt( output, output->GetLargestPossibleRegion() );
    for ( ; !outputIt.IsAtEnd() && !resultIt.IsAtEnd() ; ++outputIt, ++resultIt )
    {
      outputIt.Set( !(resultIt.Get() || outputIt.Get()) );
      progress.CompletedPixel();
    }
  }
  
 
   template <class TInputImage, class TOutputImage>
   void 
   itkLeftAtriumIsolationFilter<TInputImage, TOutputImage>
   ::PrintSelf( std::ostream &os, itk::Indent indent ) const
   {
     Superclass::PrintSelf(os, indent);
     os << indent << "FilterRadius: " << this->GetFilterRadius() << std::endl
        << indent << "No Value: " << this->GetNoValue() << std::endl
        << indent << "Interest Point Value: " << this->GetInterestPointValue() << std::endl;
   }
 
} // end namespace itk

#endif
