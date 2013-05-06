#ifndef RESAMPLEVOLUME_H
#define RESAMPleVOLUME_H

#include <stdlib.h>
#include <stdio.h>
// Include the header files for the ResampleImageFilter and the Gaussian smoothing filter.
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
// The resampling filter will need a Transform in order to map point
// coordinates and will need an interpolator in order to compute intensity
// values for the new resampled image. In this particular case we use the
// \doxygen{IdentityTransform} because the image is going to be resampled by
// preserving the physical extent of the sampled region. The Linear
// interpolator is used as a common trade-off, although arguably we should use
// one type of interpolator for the in-plane subsampling process and another
// one for the inter-slice supersampling, but again, one should wonder why to
// enter into technical sophistication here, when what we are doing is to
// cover-up for an improper acquisition of medical data, and we are just trying
// to make it look as if it was correctly acquired.
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"

using namespace std;

// Defining image types
const unsigned int Dimension = 3;
typedef itk::Image<float, Dimension> ImageType;
typedef itk::Image<int, Dimension> ImageTypeIdx;
typedef float InputPixelType;
typedef float  InternalPixelType;
typedef itk::Image< InternalPixelType, Dimension >   InternalImageType;
typedef itk::RecursiveGaussianImageFilter<ImageType, InternalImageType > GaussianFilterType;
typedef itk::ResampleImageFilter<InternalImageType, ImageType >  ResampleFilterType;
typedef vector<float> vec1f;

ImageType::Pointer ResampleVolumeToBe1Spacing(ImageType::ConstPointer IpVolume, vec1f inputSpacing, vec1f isoSpacing);

#endif // RESAMPLEVOLUME_H
