import astropy.units as u
from utils.cube_mask import cube_mask, cube_mask_channel
import numpy as np
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
u.add_enabled_units(u.def_unit(["K (Ta*)"], represents=u.K))

def mom0_normal(cube, vrange):
    vlow, vhigh = vrange
    mom0 = cube.spectral_slab(vlow*u.km/u.s,vhigh*u.km/u.s).moment0()
    return mom0
    
# in this case, pixels below threshold are masked out
def mom0_sigma_mask(cube, vrange, threshold):
    vlow, vhigh = vrange
    cube = cube.spectral_slab(vlow*u.km/u.s,vhigh*u.km/u.s)
    cube_masked = cube_mask(cube, threshold)
    mom0 = np.nansum(cube_masked, axis=0)
    return mom0
    
# in this case, channels below threshold are masked out
def mom0_sigma_channel_mask(cube, vrange, threshold):
    vlow, vhigh = vrange
    # cube = cube.spectral_slab(vlow*u.km/u.s,vhigh*u.km/u.s)
    cube_masked = cube_mask_channel(cube, threshold)
    mom0 = np.nansum(cube_masked, axis=0)
    return mom0
    
def mom1_sigma_channel_mask(cube, vrange, threshold):
    vlow, vhigh = vrange
    # cube = cube.spectral_slab(vlow*u.km/u.s,vhigh*u.km/u.s)
    cube_masked = cube_mask_channel(cube, threshold)
    v_cube = np.reshape(cube.spectral_axis.value, (-1,1,1))
    mom0 = np.nansum(cube_masked, axis=0)
    mom1 = np.nansum(cube_masked * v_cube, axis=0) / mom0
    return mom1
    
def mom2_sigma_channel_mask(cube, vrange, threshold):
    vlow, vhigh = vrange
    # cube = cube.spectral_slab(vlow*u.km/u.s,vhigh*u.km/u.s)
    cube_masked = cube_mask_channel(cube, threshold)
    v_cube = np.reshape(cube.spectral_axis.value, (-1,1,1))
    mom0 = np.nansum(cube_masked, axis=0)
    mom1 = np.nansum(cube_masked * v_cube, axis=0) / mom0
    mom2 = np.sqrt(np.nansum(cube_masked * (v_cube - mom1)**2, axis=0) / mom0)
    return mom2
    
if __name__ == "__main__":
    cube_path = input("input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
    print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")
    vlims = input("Input vmin and vmax in km/s separated by ',' to make the moment maps: ")
    vrange = list(map(float,vlims.split(',')))

    mom0_map = mom0_normal(cube,vrange)

    vmax = np.nanpercentile(mom0_map.value,99)

    fig = plt.figure()  

    axis_format = input('format the axes in pixel or world coordinates?: ')
    
    if axis_format == 'world':
        wcs = cube.wcs.celestial
        ax = fig.add_subplot(111, projection=wcs)
    elif axis_format == 'pixel':
        ax = fig.add_subplot(111)
    else:
        print('unrecognised input, default to pixel coordinates')
        ax = fig.add_subplot(111)

    im = ax.imshow(mom0_map.value, origin='lower', vmin=0, vmax=vmax)
    plt.colorbar(im, ax=ax, label=mom0_map.unit)
    title = input('customise title? enter to skip: ')
    if title == '':
        ax.set_title('Mom0 Map')
    else:
        ax.set_title(title)
    plt.show()
    save_path = input('save the figure? type the path, or enter to skip: ')
    if save_path != '':
        fig.savefig(save_path, dpi=150, bbox_inches='tight', pad_inches=0)
        print(f'figure saved to {save_path}')
    else:
        print(f'figure not saved') 
