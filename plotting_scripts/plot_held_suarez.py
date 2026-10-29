"""
Plots the Held-Suarez test case.

Plots include:
- Initial fields: lon/lat slices at surface of temperature and zonal wind,
                  and lat/z slices at lon=0 of temperature and zonal wind
- Instantaneous fields: lon/z slices of temperature and zonal wind after 100
                        and 1200 days
- Instantaneous horizontal slices: lon/lat slices of temperature and pressure
                        at the surface, and temperature and zonal wind at
                        15 km height
- Zonal averages: zonal averages of temperature and zonal wind
"""

from os.path import abspath, dirname
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_ax, extract_gusto_coords, extract_gusto_field,
    tomplot_field_title, zonal_average
)

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
test = 'held_suarez'
results_file_name = f'{abspath(dirname(__file__))}/../results/{test}/field_output.nc'
plot_stem = f'{abspath(dirname(__file__))}/../figures/{test}'

all_contours = [
    np.linspace(150, 320, 18),
    np.linspace(-40, 40, 21)
]

# 1D grids for vertical regridding
coords_lon_1d = np.linspace(-180, 180, 50)
coords_lat_1d = np.linspace(-90, 90, 50)
# Dictionary to hold plotting grids -- keys are "slice_along" values
plotting_grids = {'lat': coords_lon_1d, 'lon': coords_lat_1d}

initial_xlims = (-90, 90)
initial_ylims_km = (0, 30)
yticks = [-90, 90]
ytick_labels = ['90S', '90N']

# ---------------------------------------------------------------------------- #
# Instantaneous field plot details
# ---------------------------------------------------------------------------- #
# lon/z slices (at lat=0) of temperature and zonal wind, after 100 and 1200 days
instant_field_names = ['Temperature', 'u_zonal']
instant_titles = ['Temperature', 'Zonal Wind']
instant_field_labels = [r'$T$ (K)', r'$u$ (m/s)']
instant_colour_schemes = ['OrRd', 'RdBu_r']
instant_slice_at = 0.0

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
contour_method = 'tricontour'

set_tomplot_style(fontsize=18)

data_file = Dataset(results_file_name, 'r')

# ---------------------------------------------------------------------------- #
# Zonal averages
# ---------------------------------------------------------------------------- #
zonal_field_names = ['Temperature_average', 'u_zonal_average']
time_idxs = range(2, 13)

fig, axarray = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

time_start = data_file['time'][time_idxs[0]]
time_start_in_days = time_start / (24*60*60)
time = data_file['time'][time_idxs[-1]]
time_in_days = time / (24*60*60)

for col_idx, (field_name, title, field_label, colour_scheme, contours) in \
    enumerate(zip(zonal_field_names, instant_titles, instant_field_labels,
                  instant_colour_schemes, all_contours)):
    ax = axarray[col_idx]
    # ------------------------------------------------------------------------ #
    # Data extraction and zonal averaging
    # ------------------------------------------------------------------------ #

    # Time average fields
    field_full = extract_gusto_field(data_file, field_name, time_idxs[0])
    for time_idx in time_idxs[1:]:
        field_full += extract_gusto_field(data_file, field_name, time_idx)
    field_full /= len(time_idxs)

    coords_X_full, coords_Y_full, coords_Z_full = \
        extract_gusto_coords(data_file, field_name)

    # Bins with no data points (e.g. near the poles, for fields whose mesh
    # points don't extend all the way to the poles) are omitted by
    # zonal_average, rather than being left as NaN
    zonal_mean, lat_bin_centres, level_heights = zonal_average(
        field_full, coords_X_full, coords_Y_full, coords_Z_full, num_bins=60
    )
    # Convert height coordinate from m to km
    level_heights_km = level_heights / 1000.0
    # Zonal-mean data is structured (lat x height), so build a meshgrid and
    # use the "contour" method rather than "tricontour"
    coords_lat_2d, coords_height_2d = np.meshgrid(
        lat_bin_centres, level_heights_km
    )
    # ------------------------------------------------------------------------ #
    # Plot data
    # ------------------------------------------------------------------------ #
    cmap, lines = tomplot_cmap(contours, colour_scheme)
    cf, _ = plot_contoured_field(
        ax, coords_lat_2d, coords_height_2d, zonal_mean.T, 'contour',
        contours, cmap=cmap, line_contours=lines
    )
    add_colorbar_ax(ax, cf, field_label, location='bottom', pad=0.1,
                    cbar_labelpad=-5, cbar_format='.0f')

    tomplot_field_title(ax, title, fontsize='17.0', minmax=False)
    ax.set_xlim(initial_xlims)
    ax.set_xticks(yticks)
    ax.set_xticklabels(ytick_labels, fontsize='15.0')
    ax.set_xlabel('Latitude', labelpad=-10, fontsize='15.0')
    ax.set_ylim(initial_ylims_km)
    ax.set_yticks(initial_ylims_km)
    ax.set_yticklabels(initial_ylims_km, fontsize='15.0')
    if col_idx == 0:
        ax.set_ylabel('Height (km)', labelpad=-10, fontsize='15.0')

fig.suptitle(rf'Held-Suarez: zonal averages over {int(time_in_days) - int(time_start_in_days)} days', y=0.98, fontsize=24)
# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
plot_name = f'{plot_stem}_zonal_average.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight', dpi=300)
plt.close()
