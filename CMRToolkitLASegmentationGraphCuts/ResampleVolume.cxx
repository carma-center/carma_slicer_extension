#include "ResampleVolume.h"

ImageType::Pointer ResampleVolumeToBe1Spacing(ImageType::ConstPointer IpVolume, vec1f inputSpacing, vec1f isoSpacing){
    // We create two instances of the smoothing filter, one will smooth along the
    // $X$ direction while the other will smooth along the $Y$ direction. They are
    // connected in a cascade in the pipeline, while taking their input from the
    // intensity windowing filter. Note that you may want to skip the intensity
    // windowing scale and simply take the input directly from the reader.
    GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
    GaussianFilterType::Pointer smootherY = GaussianFilterType::New();

    smootherX->SetInput(IpVolume);
    smootherY->SetInput(smootherX->GetOutput());

    // We take the image from the input and then request its array of pixel spacing values.
    ImageType::ConstPointer inputImage = IpVolume;

    //smootherX->SetSigma( isoSpacing );
    //smootherY->SetSigma( isoSpacing );
    smootherX->SetSigma( 0.5 );
    smootherY->SetSigma( 0.5);
    // We instruct the smoothing filters to act along the $X$ and $Y$ direction
    // respectively. And define the settings for avoiding the loss of intensity as
    // a result of the diffusion process that is inherited from the use of a
    // Gaussian filter.
    smootherX->SetDirection( 0 );
    smootherY->SetDirection( 1 );

    smootherX->SetNormalizeAcrossScale( true );
    smootherY->SetNormalizeAcrossScale( true );

    // Now that we have taken care of the smoothing in-plane, we proceed to
    // instantiate the resampling filter that will reconstruct an isotropic image.
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    // The resampling filter requires that we provide a Transform, that in this
    // particular case can simply be an identity transform.
    typedef itk::IdentityTransform< double, Dimension >  TransformType;

    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();

    resampler->SetTransform( transform );

    // The filter also requires an interpolator to be passed to it. In this case we
    // chose to use a linear interpolator.
    typedef itk::LinearInterpolateImageFunction<
    InternalImageType, double >  InterpolatorType;

    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    resampler->SetInterpolator( interpolator );
    resampler->SetDefaultPixelValue( 1 ); // highlight regions without source

    // The pixel spacing of the resampled dataset is loaded in a \code{SpacingType}
    // and passed to the resampling filter.
    ImageType::SpacingType spacing;
    spacing[0] = isoSpacing[0];
    spacing[1] = isoSpacing[1];
    spacing[2] = isoSpacing[2];

    resampler->SetOutputSpacing(spacing);
    // The origin and orientation of the output image is maintained, since we
    // decided to resample the image in the same physical extent of the input
    // anisotropic image.
    resampler->SetOutputOrigin(inputImage->GetOrigin());
    resampler->SetOutputDirection(inputImage->GetDirection());
    //
    // The number of pixels to use along each dimension in the grid of the
    // resampled image is computed using the ratio between the pixel spacings of the
    // input image and those of the output image. Note that the computation of the
    // number of pixels along the $Z$ direction is slightly different with the
    // purpose of making sure that we don't attempt to compute pixels that are
    // outside of the original anisotropic dataset.
    //
    ImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();

    typedef ImageType::SizeType::SizeValueType SizeValueType;

    const double dx = (inputSize[0]) * inputSpacing[0] / isoSpacing[0];
    const double dy = (inputSize[1]) * inputSpacing[1] / isoSpacing[1];
    const double dz = (inputSize[2]-1) * inputSpacing[2] / isoSpacing[2];
    // Finally the values are stored in a \code{SizeType} and passed to the the
    // resampling filter. Note that this process requires a casting since the
    // computation are performed in \code{double}, while the elements of the
    // \code{SizeType} are integers.
    ImageType::SizeType size;

    size[0] = static_cast<SizeValueType>( dx );
    size[1] = static_cast<SizeValueType>( dy );
    if(isoSpacing[2] == 1)
        size[2] = 107;
    else
        size[2] = 44;
//    size[2] = static_cast<SizeValueType>( dz );

    resampler->SetSize( size );

    // Our last action is to take the input for the resampling image filter from
    // the output of the cascade of smoothing filters, and then to trigger the
    // execution of the pipeline by invoking the \code{Update()} method on the
    // resampling filter.
    resampler->SetInput( smootherY->GetOutput() );
    resampler->Update();

    // Writing out resampled image to image pointer
    ImageType::Pointer ResampledVolume = resampler->GetOutput();
    return ResampledVolume;
}
