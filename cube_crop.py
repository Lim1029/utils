from spectral_cube import SpectralCube
import astropy.units as u
from regions import RectangleSkyRegion
from astropy.coordinates import SkyCoord
import os

import argparse
parser = argparse.ArgumentParser(description="Cube Crop Spatial")
parser.add_argument("input_path", type=str, help="Input FITS Cube Path")
parser.add_argument("output_path", type=str, help="Output FITS Cube Path")
parser.add_argument("width", type=float, help="Width in arcmin")
parser.add_argument("height", type=float, help="Height in arcmin")
parser.add_argument("--center_l", type=float, default=28.79, help="Center Longitude in degree")
parser.add_argument("--center_b", type=float, default=3.52, help="Center Latitude in degree")

args = parser.parse_args()

os.makedirs(os.path.dirname(args.output_path), exist_ok=True)


input_cube = SpectralCube.read(args.input_path)
output_path = args.output_path

center = SkyCoord(args.center_l, args.center_b, unit='deg', frame='galactic')
width = args.width * u.arcmin
height = args.height * u.arcmin

region = RectangleSkyRegion(center=center, width=width, height=height)
input_subcube = input_cube.subcube_from_regions([region])
input_subcube.write(output_path, overwrite=True)