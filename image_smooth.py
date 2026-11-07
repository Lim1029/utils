import numpy as np
from astropy.convolution import convolve, Gaussian2DKernel
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from reproject import reproject_interp
import matplotlib.pyplot as plt

def smooth_image(data, target_reso, current_reso, pixel_size):
    kernel_fwhm = np.sqrt(target_reso**2 - current_reso**2)
    kernel_stddev = kernel_fwhm / (2*np.sqrt(2*np.log(2)))
    kernel_stddev_pixel = kernel_stddev / pixel_size
    # print(kernel_stddev_pixel)
    gauss = Gaussian2DKernel(x_stddev=kernel_stddev_pixel, y_stddev=kernel_stddev_pixel)
    result = convolve(data, gauss)
    return result
    
def reproject_image(data, header, fov, target_beam):

    pixsize = (target_beam*u.arcsec / 3).to(u.deg)  # 1/3 of the beam in degrees
    naxis = round((fov*u.arcmin).to(u.deg)/pixsize)

    # configure target header 
    tgthdr = header.copy()
    tgthdr["CDELT1"] = -pixsize.value  # negative in longitude
    tgthdr["CDELT2"] = pixsize.value
    tgthdr["NAXIS1"] = int(naxis.value)
    tgthdr["NAXIS2"] = int(naxis.value)
    tgthdr["CRPIX1"] = (naxis.value+1)/2 # reference x pixel
    tgthdr["CRPIX2"] = (naxis.value+1)/2 # reference y pixel
    tgthdr['CRVAL1'] = 28.79 # b value at crpix1, hardcoded unfortunately
    tgthdr['CRVAL2'] = 3.52 # l value at crpix2, hardcoded unfortunately

    output_wcs = WCS(tgthdr)
    shape_out = (int(naxis.value),int(naxis.value))

   # because reproject_interp only accepts HDUList or HDU as input, we need to convert the Cutout2D object to HDU object
    # print(header)
    hdu = fits.PrimaryHDU(data, header=header)
    array, footprint = reproject_interp(input_data=hdu, 
                                        output_projection=output_wcs, 
                                        shape_out=shape_out)
    return array,output_wcs
    
if __name__ == "__main__":
  
    fits_path = input("input path to fits image: ")
    hdu = fits.open(fits_path)[0]    
    target_reso = float(input("input target reso (arcsec): "))
    current_reso = float(input("input current reso (arcsec): "))
    current_pixelsize = hdu.header['CDELT2']*3600 # deg2arcsec
    smoothed = smooth_image(hdu.data, target_reso, current_reso, current_pixelsize)
    plt.imshow(smoothed, origin='lower')
    plt.show()
    fov = float(input("input fov (arcmin): "))
    array,output_wcs = reproject_image(smoothed, hdu.header, fov, target_reso)
    
    # create a new hdu and hdul and export to FITS
    new_hdu = fits.PrimaryHDU(data=array, header=output_wcs.to_header())
    new_hdul = fits.HDUList([new_hdu])
    output_path = input("input output path: ")
    new_hdul.writeto(output_path, overwrite=True)