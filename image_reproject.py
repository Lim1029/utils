from astropy.io import fits
from reproject import reproject_interp


def image_reproject(data, wcs, output_wcs, shape_out):
   # because reproject_interp only accepts HDUList or HDU as input, we need to convert the Cutout2D object to HDU object
    hdu = fits.PrimaryHDU(data, header=wcs.to_header())
    array, footprint = reproject_interp(input_data=hdu, 
                                        output_projection=output_wcs, 
                                        shape_out=shape_out)
    return array