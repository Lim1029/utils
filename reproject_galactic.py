############### User defined #################### 
starlink_dir = '/home/mingkang/star-2025A'
############### User defined #################### 

from astropy.io import fits
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.wcs import WCS

import argparse
parser = argparse.ArgumentParser(description="Image Reprojection")
parser.add_argument("input_path", type=str, help="Input FITS Image Path")
parser.add_argument("output_path", type=str, help="Output FITS Image Path")
args = parser.parse_args()

file_path = args.input_path #input fits to be reprojected
output_path = args.output_path #output fits name

output_frame = 'galactic' # or icrs

# Do not modify anything below

wcs = WCS(fits.getheader(file_path)).celestial

# Use find_optimal_celestial_wcs with the 2D data
wcs_out, shape_out = find_optimal_celestial_wcs(wcs, frame=output_frame)

import numpy as np
# Create a reference FITS file with the WCS and shape
hdu = fits.PrimaryHDU(data=np.zeros(shape_out), header=wcs_out.to_header())
hdul = fits.HDUList([hdu])

# Save the reference FITS file
reference_fits_path = './reference_wcs.fits'
hdul.writeto(reference_fits_path, overwrite=True)

print(f"reference FITS file created at: {reference_fits_path}")

import subprocess

subprocess.run(['bash', '-c', f'''
export STARLINK_DIR={starlink_dir}
source $STARLINK_DIR/etc/profile
convert
kappa
wcsalign {file_path} {output_path} ref={reference_fits_path} accept
rm {reference_fits_path}
'''])

