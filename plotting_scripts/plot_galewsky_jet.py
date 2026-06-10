"""
Plots the relative vorticity field from the Galewsky jet test case, comparing
an icosahedral sphere and cubed sphere mesh, at 6 days.
"""
from os.path import abspath, dirname
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field, add_colorbar_fig,
    tomplot_field_title, extract_gusto_coords, extract_gusto_field,
    plot_icosahedral_sphere_panels, plot_cubed_sphere_panels,
    regrid_horizontal_slice
)

results_dirs = ['galewsky_jet_cubedsphere', 'galewsky_jet_icosahedralsphere']

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
results_stem = f'{abspath(dirname(__file__))}/../results/'
plot_stem = f'{abspath(dirname(__file__))}/../figures/galewsky_jet'

# ---------------------------------------------------------------------------- #
# Plot details
# ---------------------------------------------------------------------------- #
subplots_x = 2
subplots_y = 1
field_name = 'RelativeVorticity'
time_idx = -1
subtitles = ['Cubed-Sphere', 'Icosahedral Sphere']

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #

spherical_centre = (-15, 75)
projection = ccrs.Orthographic(spherical_centre[0], spherical_centre[1])
contours = np.linspace(-2e-4, 2e-4, 17)
colour_scheme = 'RdBu_r'
field_label = r'$\zeta$ (s$^{-1}$)'
contour_method = 'contour'  # Best method for orthographic projection
# We need to regrid onto lon-lat grid -- specify that here
lon_1d = np.linspace(-180.0, 180.0, 500)
lat_1d = np.linspace(-90, 90, 500)
coords_lon, coords_lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')

# Things that are likely the same for all plots --------------------------------
set_tomplot_style(20)

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #

fig = plt.figure(figsize=(15, 6))

# Save filled contour data for adding colourbars at end
all_cf = []

for i, (subtitle, results_dir) in enumerate(zip(subtitles, results_dirs)):

    results_file_name = f'{results_stem}{results_dir}/field_output.nc'
    data_file = Dataset(results_file_name, 'r')

    ax = fig.add_subplot(subplots_y, subplots_x, 1+i, projection=projection)

    # Data extraction ----------------------------------------------------------
    orig_field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
    orig_coords_lon, orig_coords_lat = extract_gusto_coords(data_file, field_name)
    time = data_file['time'][time_idx] / (24.*60.*60.)

    field_data = regrid_horizontal_slice(coords_lon, coords_lat, orig_coords_lon,
                                         orig_coords_lat, orig_field_data,
                                         periodic_fix='sphere')

    cmap, _ = tomplot_cmap(contours, colour_scheme)
    _, lines = tomplot_cmap(contours, colour_scheme, remove_contour=0.0)

    # Plot data ----------------------------------------------------------------

    cf, _ = plot_contoured_field(
        ax, coords_lon, coords_lat, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, projection=projection
    )
    all_cf.append(cf)

    tomplot_field_title(
        ax, subtitle, minmax=False, field_data=field_data
    )

    if i == 0:
        plot_cubed_sphere_panels(ax)
    else:
        plot_icosahedral_sphere_panels(ax)

# Colour bars need adding at end when the full figure exists -------------------
fig.subplots_adjust(wspace=0.12, hspace=-0.05)

cbar_labelpad = -70
data_format = '.0e'

cb = add_colorbar_fig(
    fig, all_cf[i], field_label, ax_idxs=[1], location='right',
    cbar_labelpad=cbar_labelpad, cbar_format=data_format
)

# Change cbar tick labels to scientific notation
cbar_ticklabels = [r'$-2\times 10^{-4}$', r' $2\times 10^{-4}$']
cb.ax.set_yticklabels(cbar_ticklabels)

# Save figure ------------------------------------------------------------------
fig.suptitle(rf'Unstable Jet, $t=$ {int(time)} days', y=0.98, fontsize=20)
plot_name = f'{plot_stem}.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight', dpi=300)
plt.close()
