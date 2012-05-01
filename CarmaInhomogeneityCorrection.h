
#ifndef InhomogeneityCorrectionFilter_H
#define InhomogeneityCorrectionFilter_H

// std includes
#include <algorithm>
#include <cmath>
#include <limits>

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
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkCommand.h>
#include <itkLinearInterpolateImageFunction.h>
 

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>

  
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
  class InhomogeneityCorrectionFilter :
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
      typedef InhomogeneityCorrectionFilter Self;
      typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
      typedef SmartPointer<Self> Pointer;
      typedef SmartPointer<const Self> ConstPointer;

      /** Method for creation through the object factory. */
      itkNewMacro(Self);

      /** Run-time type information (and related methods). */
      itkTypeMacro(InhomogeneityCorrectionFilter, ImageToImageFilter);
  
      /** Image typedef support. */
      typedef typename InputImageType::PixelType  InputPixelType;
      typedef typename OutputImageType::PixelType OutputPixelType;
    
      typedef typename InputImageType::RegionType  InputImageRegionType;
      typedef typename OutputImageType::RegionType OutputImageRegionType;

      typedef typename InputImageType::SizeType InputSizeType;

      typedef typename InputImageType::PointType PointType;

      typedef typename itk::LinearInterpolateImageFunction<InputImageType> LinearInterpolatorType;


      typedef itk::ImageRegionConstIteratorWithIndex<InputImageType> ConstIteratorType;


      //// get/set-ers:
      itkSetMacro(SampleSpacing, float); // millimeters
      itkGetConstReferenceMacro(SampleSpacing, float); // millimeters

      itkSetMacro(PolynomialOrder, int); // millimeters
      itkGetConstReferenceMacro(PolynomialOrder, int); // millimeters

     /** InhomogeneityCorrectionFilter needs a larger input requested region than
       * the output requested region.  As such, InhomogeneityCorrectionFilter needs
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
    InhomogeneityCorrectionFilter();
    virtual ~InhomogeneityCorrectionFilter() {}

    void PrintSelf(std::ostream &os, itk::Indent indent) const;

    void GenerateData();

  private:
    InhomogeneityCorrectionFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector<std::pair<PointType, double> > m_Samples;

    float m_SampleSpacing;
    int m_PolynomialOrder;

    typename LinearInterpolatorType::Pointer m_DataInterpolator;
    typename LinearInterpolatorType::Pointer m_MaskInterpolator;


    void SampleImage(itk::ProgressReporter *progress);


    vnl_matrix<double> Eval2dPoly(vnl_matrix<double> x, vnl_matrix<double> y, vnl_matrix<double> z, vnl_matrix<double> coeffs);
   
 }; // end class InhomogeneityCorrectionFilter


  // DEFINITIONS:

  template <class TInputImage, class TOutputImage>
  InhomogeneityCorrectionFilter<TInputImage, TOutputImage>
  ::InhomogeneityCorrectionFilter()
  :m_SampleSpacing(1.0f),
  m_PolynomialOrder(2)
  {
    // constructor
  }

  template <class TInputImage, class TOutputImage>
  vnl_matrix<double> InhomogeneityCorrectionFilter<TInputImage, TOutputImage>
    ::Eval2dPoly( vnl_matrix<double> x, vnl_matrix<double> y, vnl_matrix<double> z, vnl_matrix<double> coeffs )
  {
/*
   [sizexR, sizexC] = size(x);
   [sizeyR, sizeyC] = size(y);
   [sizezR, sizezC] = size(z);
   numVals = sizexR;

   order = 0.5 * (sqrt(8*length(coeffs)+1) - 5);

   zbar = zeros(numVals,1);
   column = 1;
   for xpower = 0:order
     for ypower = 0:order
       for zpower = 0:order
         if (xpower + ypower + zpower <= order)
           zbar = zbar + (coeffs(column) .* x.^xpower .* y.^ypower.*z.^zpower);
           column = column + 1;
         end
       end
     end
   end
*/

    int sizexR = x.rows();
    int sizexC = x.cols();
    int sizeyR = y.rows();
    int sizeyC = y.cols();
    int sizezR = z.rows();
    int sizezC = z.cols();
    int numVals = sizeyR;

    int size_coeffs = coeffs.size();
    double order = 0.5 * (std::sqrt((double)(8*size_coeffs+1)) - 5);
  
    vnl_matrix<double> zbar(numVals, 1, 0.0);
    int column = 0;
    for (double xpower = 0; xpower <= m_PolynomialOrder; xpower++) {
      for (double ypower = 0; ypower <= m_PolynomialOrder; ypower++) {
        for (double zpower = 0; zpower <= m_PolynomialOrder; zpower++) {
          if (xpower + ypower + zpower <= m_PolynomialOrder) {
            for (int i=0; i < numVals; i++) {
              zbar(i,0) = zbar(i,0) + (coeffs(column,0) * pow(x(i,0),xpower) * pow(y(i,0),ypower) * pow(z(i,0),zpower));
            }  
            column++;
          }     
        }
      }
    }
    return zbar;
  }


 
  template< class TInputImage, class TOutputImage>
  void
  InhomogeneityCorrectionFilter< TInputImage, TOutputImage>
  ::GenerateData()
  {
    // Allocate output
    typename OutputImageType::Pointer output = this->GetOutput();
    output->SetBufferedRegion( output->GetRequestedRegion() );   
    output->Allocate();
    output->FillBuffer ( itk::NumericTraits<OutputPixelType>::Zero ); 

    // Get Input pointers
    typename InputImageType::ConstPointer input_data = this->GetInput(0);
    typename InputImageType::ConstPointer input_endo = this->GetInput(1);
    typename InputImageType::ConstPointer input_roi = this->GetInput(2);


    if (InputImageDimension != 3) {
      itkExceptionMacro ("InhomogeneityCorrectionFilter only works in 3D");
    }

    // support progress methods/callbacks
    unsigned int progress_sections = 16;
    itk::ProgressReporter progress(this, 0, progress_sections);


    typename InputImageType::PointType origin = input_data->GetOrigin();
    typename InputImageType::RegionType region = input_data->GetLargestPossibleRegion();
    typename InputImageType::SpacingType spacing = input_data->GetSpacing();
    typename InputImageType::SizeType size = region.GetSize();

    m_DataInterpolator = LinearInterpolatorType::New();
    m_DataInterpolator->SetInputImage(input_data);

    m_MaskInterpolator = LinearInterpolatorType::New();
    m_MaskInterpolator->SetInputImage(input_endo);

    SampleImage(&progress);
    progress.CompletedPixel();

    int num_points = m_Samples.size();
    //std::cout << "Found " << m_Samples.size() << " points\n";

    if (num_points < 1) {
      std::cerr << "InhomogeneityCorrectionFilter: Error: did not find any sampled points!\n";
      itkExceptionMacro ("Error, did not find any sampled points! (reduce sampling rate?)");
    }


    //number of combinations of coefficients in resulting polynomial
    int numCoeffs = (m_PolynomialOrder+3)*(m_PolynomialOrder+2)/2;

    vnl_matrix<double> pts(num_points, numCoeffs, 0.0);

    vnl_matrix<double> x_array(num_points, 1, 0.0);
    vnl_matrix<double> y_array(num_points, 1, 0.0);
    vnl_matrix<double> z_array(num_points, 1, 0.0);
    vnl_matrix<double> a_array(num_points, 1, 0.0);

    for (int i=0; i < num_points; i++) {
      PointType p = m_Samples[i].first;
      x_array(i,0) = p[0];
      y_array(i,0) = p[1];
      z_array(i,0) = p[2];
      a_array(i,0) = m_Samples[i].second;
    }

    double scale_x = 1.0 / x_array.absolute_value_max();
    double scale_y = 1.0 / y_array.absolute_value_max();
    double scale_z = 1.0 / z_array.absolute_value_max();
    double scale_a = 1.0 / a_array.absolute_value_max();

    vnl_matrix<double> xs = x_array * scale_x;
    vnl_matrix<double> ys = y_array * scale_y;
    vnl_matrix<double> zs = z_array * scale_z;
    vnl_matrix<double> as = a_array * scale_a;

    int column = 0;
    for (double xpower = 0; xpower <= m_PolynomialOrder; xpower++) {
      for (double ypower = 0; ypower <= m_PolynomialOrder; ypower++) {
        for (double zpower = 0; zpower <= m_PolynomialOrder; zpower++) {
          if (xpower + ypower + zpower <= m_PolynomialOrder) {
            for (int i=0; i < num_points; i++) {
              pts(i, column) = std::pow(xs(i,0),xpower) * std::pow(ys(i,0),ypower) * std::pow(zs(i,0),zpower);
            }
            column++;
          }       
        }
      }
    }

    vnl_svd<double> svd(pts);
    progress.CompletedPixel();

    vnl_matrix<double> u = svd.U();
    vnl_matrix<double> s = svd.W();
    vnl_matrix<double> v = svd.V();

    //std::cerr << v << "\n";

    double eps = std::numeric_limits<double>::epsilon();
    double sigma = std::pow (eps, 1.0/(double)m_PolynomialOrder);

    vnl_matrix<double> qqs(num_points, numCoeffs, 0.0);

    for (int i=0; i < numCoeffs; i++) {
      if (s(i,i) >= sigma) {
        qqs(i,i) = 1.0 / s(i,i);
      } else {
        qqs(i,i) = 0.0;
      }
    }

    //vnl_matrix<double> coeffs = v * qqs.transpose() * u.transpose() * as;
    vnl_matrix<double> coeffs = svd.solve(as);

    progress.CompletedPixel();

    //  % scale the coefficients so they are correct for the unscaled data
    //  column = 1;
    //  for xpower = 0:order
    //    for ypower = 0:order
    //      for zpower = 0:order
    //        if (xpower + ypower + zpower <= order)
    //          coeffs(column) = coeffs(column) * scalex^xpower * scaley^ypower * scalez^zpower / scalea;
    //          column = column + 1;
    //        end
    //      end
    //    end
    //  end

    column = 0;
    for (double xpower = 0; xpower <= m_PolynomialOrder; xpower++) {
      for (double ypower = 0; ypower <= m_PolynomialOrder; ypower++) {
        for (double zpower = 0; zpower <= m_PolynomialOrder; zpower++) {
          if (xpower + ypower + zpower <= m_PolynomialOrder) {
            coeffs(column,0) = coeffs(column,0) * pow(scale_x,xpower) * pow(scale_y,ypower) * pow(scale_z, zpower) / scale_a;
            column++;
          }       
        }
      }
    }

    // coeffs = coeffs./mean(ii(roiBP > 0));
    coeffs = coeffs / a_array.mean();

    // collect together roi points
    ConstIteratorType data_it(input_data, input_data->GetLargestPossibleRegion());
    ConstIteratorType roi_it(input_roi, input_roi->GetLargestPossibleRegion());
    data_it.GoToBegin();
    roi_it.GoToBegin();
    std::vector<std::pair<PointType, double> > roi_points;
    for ( ; !data_it.IsAtEnd(); ++data_it,++roi_it) {
      typename InputImageType::PixelType image_value = data_it.Get();
      typename InputImageType::PixelType roi_value = roi_it.Get();
      if (roi_value > 0) {
        typename InputImageType::IndexType index = roi_it.GetIndex();
        PointType p;
        input_roi->TransformIndexToPhysicalPoint(index,p);
        roi_points.push_back(std::make_pair(p,image_value));
      }
    }

    vnl_matrix<double> roi_x(roi_points.size(), 1, 0.0);
    vnl_matrix<double> roi_y(roi_points.size(), 1, 0.0);
    vnl_matrix<double> roi_z(roi_points.size(), 1, 0.0);
    vnl_matrix<double> roi_a(roi_points.size(), 1, 0.0);

    for (int i=0; i < roi_points.size(); i++) {
      PointType p = roi_points[i].first;
      roi_x(i,0) = p[0];
      roi_y(i,0) = p[1];
      roi_z(i,0) = p[2];
      roi_a(i,0) = roi_points[i].second;
    }

    progress.CompletedPixel();

    vnl_matrix<double> correction = Eval2dPoly(roi_x, roi_y, roi_z, coeffs);

    progress.CompletedPixel();


    itk::ImageRegionIterator< OutputImageType > output_it( output , output->GetLargestPossibleRegion() );

    int index = 0;
    data_it.GoToBegin();
    roi_it.GoToBegin();
    for ( ; !output_it.IsAtEnd(); ++output_it, ++roi_it, ++data_it )
    {
      typename InputImageType::PixelType image_value = data_it.Get();
      typename InputImageType::PixelType roi_value = roi_it.Get();

      if (roi_value > 0) {
        double correction_value = correction(index,0);
        if (correction_value < 0.3) {
          image_value = 0;
        } else {
          image_value /= correction_value;
        }

        //image_value /= correction(index, 0);
        index++;
      }

      output_it.Set(image_value);
    }

    progress.CompletedPixel();

}
  


template <class TInputImage, class TOutputImage>
void InhomogeneityCorrectionFilter<TInputImage, TOutputImage>
::SampleImage( itk::ProgressReporter *progress )
{
  typename InputImageType::ConstPointer input_data = this->GetInput(0);
  typename InputImageType::ConstPointer input_endo = this->GetInput(1);
  typename InputImageType::ConstPointer input_roi = this->GetInput(2);

  typename InputImageType::PointType origin = input_data->GetOrigin();
  typename InputImageType::RegionType region = input_data->GetLargestPossibleRegion();
  typename InputImageType::SpacingType spacing = input_data->GetSpacing();
  typename InputImageType::SizeType size = region.GetSize();


  //float margin = m_SampleSpacing * 2;

  float start_x = origin[0];
  float start_y = origin[1];
  float start_z = origin[2];
  float end_x = origin[0] + (spacing[0] * size[0]);
  float end_y = origin[1] + (spacing[1] * size[1]);
  float end_z = origin[2] + (spacing[2] * size[2]);

  float possible_samples = (end_x - start_x) / m_SampleSpacing;
  possible_samples *= (end_y - start_y) / m_SampleSpacing;
  possible_samples *= (end_z - start_z) / m_SampleSpacing;

  int sample_count = 0;

  for (float x_position = start_x; x_position < end_x; x_position += m_SampleSpacing) {
    for (float y_position = start_y; y_position < end_y; y_position += m_SampleSpacing) {
      for (float z_position = start_z; z_position < end_z; z_position += m_SampleSpacing) {

        sample_count++;

        if (sample_count > (possible_samples / 10)) {
          sample_count = 0;
          progress->CompletedPixel();
        }

        PointType p;
        p[0] = x_position;
        p[1] = y_position;
        p[2] = z_position;

        typename LinearInterpolatorType::ContinuousIndexType index;
        input_data->TransformPhysicalPointToContinuousIndex(p, index);

        typename InputImageType::PixelType mask_value = m_MaskInterpolator->EvaluateAtContinuousIndex(index);

        if (mask_value >= 1.0) {
          typename InputImageType::PixelType image_value = m_DataInterpolator->EvaluateAtContinuousIndex(index);
          m_Samples.push_back(std::make_pair(p,image_value));
        }
      }
    }
  }
}

 
template <class TInputImage, class TOutputImage>
void 
InhomogeneityCorrectionFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "SampleSpacing: " << this->GetSampleSpacing() << std::endl;
}

#endif
