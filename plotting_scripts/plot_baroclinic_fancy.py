"""
Plots the Baroclinic wave test case.

This plots:
(a) Surface Pressure , (b) Temperature
"""
from os.path import abspath, dirname
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
import cartopy.crs as ccrs
import numpy as np
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_ax, tomplot_field_title, tomplot_contours,
    extract_gusto_coords, extract_gusto_field, reshape_gusto_data,
    regrid_horizontal_slice
)

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these paths need editing, which will usually involve
# removing the abspath part to set directory paths relative to this file
results_dir = f'/media/thomas/1670402270400B47/results/baroclinic_wave'
plot_stem = f'{abspath(dirname(__file__))}/../../source/figures/'

# ---------------------------------------------------------------------------- #
# Final plot details
# ---------------------------------------------------------------------------- #
field_names = ['Pressure_Vt', 'Temperature']
slice_along = 'z'
colours = ['PRGn', 'Reds']
centre_titles = ['Surface pressure',  'Temperature']
field_labels = [r'$P \ / $hPa', r'$T \ / $K']

# We need to regrid onto lon-lat grid -- specify that here
lon_1d = np.linspace(-180.0, 180.0, 80)
lat_1d = np.linspace(-90, 90, 80)
coords_lon, coords_lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')

spherical_centre = (120, 60)
projection = ccrs.Orthographic(spherical_centre[0], spherical_centre[1])

contours = [
    np.linspace(94000, 102000, 17),
    np.linspace(230, 330, 21)
]
remove_contours = [100000.0, None]
cmap_rescales = ['top', None]

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
contour_method = 'contour'
domain_limit = {'X' : (-180, 180), 'Y' : (-90, 90)}

time_idx = 10*4
# Things that are likely the same for all plots --------------------------------
set_tomplot_style()
level = 0

pressure_ticks = [94000, 96000, 98000, 100000, 102000]
pressure_labels = ['940', '960', '980', '1000', '1020']
temp_ticks = [240, 280, 320]
temp_labels = ['240', '280', '320']

cbar_ticks = [pressure_ticks, temp_ticks]
cbar_ticklabels = [pressure_labels, temp_labels]

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #

results_file_name = f'{results_dir}/field_output.nc'
data_file = Dataset(results_file_name, 'r')
time = data_file['time'][time_idx]
time_in_days = time / (24*60*60)

fig = plt.figure(figsize=(15, 8))

for i, (field_name, field_label, contour, centre_title, colour,
        cmap_rescale, cbar_tick, cbar_ticklabel, remove_contour) in enumerate(zip(
        field_names, field_labels, contours, centre_titles,
        colours, cmap_rescales, cbar_ticks, cbar_ticklabels, remove_contours)):

    ax = fig.add_subplot(1, 2, 1+i, projection=projection)

    # Data extraction ------------------------------------------------------
    field_full = extract_gusto_field(data_file, field_name, time_idx)
    coords_X_full, coords_Y_full, coords_Z_full = extract_gusto_coords(data_file, field_name)

    # Reshape
    field_full, coords_X_full, coords_Y_full, _ = \
        reshape_gusto_data(field_full, coords_X_full, coords_Y_full, coords_Z_full)

    orig_field_data = field_full[:, level]
    orig_coords_lon = coords_X_full[:, level]
    orig_coords_lat = coords_Y_full[:, level]

    field_data = regrid_horizontal_slice(coords_lon, coords_lat, orig_coords_lon,
                                         orig_coords_lat, orig_field_data,
                                         periodic_fix='sphere')

    cmap, lines = tomplot_cmap(
        contour, colour, remove_contour=remove_contour,
        cmap_rescale_type=cmap_rescale
    )

    # Plot data ------------------------------------------------------------
    cf, _ = plot_contoured_field(
        ax, coords_lon, coords_lat, field_data, contour_method, contour,
        cmap=cmap, line_contours=lines, projection=projection
    )

    # cf, _ = plot_contoured_field(
    #     ax, coords_lon, coords_lat, field_data, contour_method, contour,
    #     cmap=cmap, plot_contour_lines=False, projection=projection
    # )

    cb = add_colorbar_ax(ax, cf, field_label, location='bottom', cbar_ticks=cbar_tick, pad=0.07)
    cb.ax.set_xticklabels(cbar_ticklabel)

    # Labels ---------------------------------------------------------------
    # Adding titles
    ax.set_title(centre_title, fontsize='17.0')

# Save figure ------------------------------------------------------------------
fig.suptitle(rf'Baroclinic wave, Day: {int(time_in_days):02d}', y=0.98, fontsize=20)
plot_name = 'baroclinic_wave_fancy.png'

print(f'Saving figure to {plot_stem}/{plot_name}')
png_name = f"{plot_stem}/{plot_name}"
fig.savefig(png_name, bbox_inches='tight')
