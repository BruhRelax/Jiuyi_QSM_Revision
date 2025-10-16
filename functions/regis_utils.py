import SimpleITK as sitk

def check_alignment(reference_image: sitk.Image, moving_image: sitk.Image) -> bool:
    """
    Checks if two images are aligned by comparing their metadata.
    Returns True if aligned, False otherwise.
    """
    is_aligned = all([
        reference_image.GetSize() == moving_image.GetSize(),
        abs(reference_image.GetOrigin()[0] - moving_image.GetOrigin()[0]) < 1e-6,
        abs(reference_image.GetSpacing()[0] - moving_image.GetSpacing()[0]) < 1e-6,
        abs(reference_image.GetDirection()[0] - moving_image.GetDirection()[0]) < 1e-6
    ])
    return is_aligned

def register_image(reference_image: sitk.Image, moving_image: sitk.Image, transform_type: str, log_callback=None):
    """
    Performs registration and returns both the registered image and the transform.

    Args:
        log_callback (function, optional): A function to log messages to the GUI.
    """
    if log_callback:
        log_callback(f"Performing {transform_type} registration...")

    if transform_type.lower() == 'rigid':
        transform = sitk.Euler3DTransform()
    elif transform_type.lower() == 'affine':
        transform = sitk.AffineTransform(reference_image.GetDimension())
    else:
        raise ValueError(f"Unknown transform type: '{transform_type}'.")

    initial_transform = sitk.CenteredTransformInitializer(
        reference_image, moving_image, transform, sitk.CenteredTransformInitializerFilter.GEOMETRY
    )
    
    registration_method = sitk.ImageRegistrationMethod()
    
    # Configure the metric
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)

    # Configure the interpolator
    registration_method.SetInterpolator(sitk.sitkLinear)
    
    # Configure the optimizer
    registration_method.SetOptimizerAsGradientDescent(
        learningRate=1.0, numberOfIterations=100,
        convergenceMinimumValue=1e-6, convergenceWindowSize=10
    )
    registration_method.SetOptimizerScalesFromPhysicalShift()
    
    # Configure the transform
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    final_transform = registration_method.Execute(
        sitk.Cast(reference_image, sitk.sitkFloat64),
        sitk.Cast(moving_image, sitk.sitkFloat64)
    )

    if log_callback:
        log_callback(f"Optimizer stop condition: {registration_method.GetOptimizerStopConditionDescription()}")
        log_callback(f"Final metric value: {registration_method.GetMetricValue():.4f}")

    registered_image = sitk.Resample(
        moving_image, reference_image, final_transform,
        sitk.sitkLinear, 0.0, moving_image.GetPixelID()
    )
    
    return registered_image, final_transform


