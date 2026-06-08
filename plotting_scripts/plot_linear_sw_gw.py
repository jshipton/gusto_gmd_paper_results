"""
Plots the linear gravity wave test case.

The initial conditions are plotted:
(a) the velocity field, (b) the depth field

And after 5 days, this plots:
(a) depth
"""
from os.path import abspath, dirname
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_ax, plot_field_quivers, tomplot_field_title,
    extract_gusto_coords, extract_gusto_field, regrid_horizontal_slice
)

test_name = 'linear_sw_gw_ssprk3_RTCF2_n36_dt15'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
# When copying this example these paths need editing, which will usually involve
# removing the abspath part to set directory paths relative to this file
results_file_name = f'results/{test_name}/field_output.nc'
plot_stem = f'results/{test_name}'

# ---------------------------------------------------------------------------- #
# Initial plot details
# ---------------------------------------------------------------------------- #
init_field_names = ['u', 'D_minus_Dbalanced']
init_colour_schemes = ['Oranges', 'YlGnBu']
init_field_labels = [r"$|u'|$ (m s$^{-1}$)", r"$D'$ (m)"]
init_contours_to_remove = [None, None, None]
init_contours = [
    np.linspace(0, 40, 21),
    np.linspace(-100, 100, 21)
]

# ---------------------------------------------------------------------------- #
# Final plot details
# ---------------------------------------------------------------------------- #
final_field_names = ["D_minus_Dbalanced"]
final_colour_schemes = ['YlGnBu']
final_field_labels = [r"$D'$ (m)"]
final_contours = [
    np.linspace(-100, 100, 21),
]
final_contours_to_remove = [None]

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
projection = ccrs.Robinson()
contour_method = 'contour'
xlims = [-180, 180]
ylims = [-90, 90]

# We need to regrid onto lon-lat grid -- specify that here
lon_1d = np.linspace(-180.0, 180.0, 120)
lat_1d = np.linspace(-90, 90, 120)
lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d, indexing='ij')

# Things that are likely the same for all plots --------------------------------
set_tomplot_style()
data_file = Dataset(results_file_name, 'r')

# ---------------------------------------------------------------------------- #
# INITIAL PLOTTING
# ---------------------------------------------------------------------------- #
fig = plt.figure(figsize=(15, 5))
time_idx = 0

for i, (field_name, colour_scheme, field_label, contour_to_remove, contours) in \
    enumerate(zip(
        init_field_names, init_colour_schemes,
        init_field_labels, init_contours_to_remove, init_contours)):

    # Make axes
    ax = fig.add_subplot(1, 2, 1+i, projection=projection)

    # Data extraction ----------------------------------------------------------
    if field_name == 'u':
        zonal_data = extract_gusto_field(data_file, 'u_zonal', time_idx=time_idx)
        meridional_data = extract_gusto_field(data_file, 'u_meridional', time_idx=time_idx)
        field_data = np.sqrt(zonal_data**2 + meridional_data**2)
        coords_X, coords_Y = extract_gusto_coords(data_file, 'u_zonal')

    else:
        field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
        coords_X, coords_Y = extract_gusto_coords(data_file, field_name)
    time = data_file['time'][time_idx] / (60.*60.)

    # Regrid onto lon-lat grid
    field_data = regrid_horizontal_slice(
        lon_2d, lat_2d, coords_X, coords_Y, field_data, periodic_fix='sphere'
    )

    # Plot data ----------------------------------------------------------------
    cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=contour_to_remove)
    cf, _ = plot_contoured_field(
        ax, lon_2d, lat_2d, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, projection=projection
    )

    add_colorbar_ax(ax, cf, field_label, location='bottom', cbar_labelpad=-10)
    tomplot_field_title(ax, None, minmax=True, field_data=field_data)

    # Add quivers --------------------------------------------------------------
    if field_name == 'u':
        # Need to re-grid to lat-lon grid to get sensible looking quivers
        regrid_zonal_data = regrid_horizontal_slice(
            lon_2d, lat_2d, coords_X, coords_Y, zonal_data,
            periodic_fix='sphere'
        )
        regrid_meridional_data = regrid_horizontal_slice(
            lon_2d, lat_2d, coords_X, coords_Y, meridional_data,
            periodic_fix='sphere'
        )
        plot_field_quivers(
            ax, lon_2d, lat_2d, regrid_zonal_data, regrid_meridional_data,
            spatial_filter_step=6, magnitude_filter=1.0,
            projection=ccrs.PlateCarree()
        )

# Save figure ------------------------------------------------------------------
fig.subplots_adjust(wspace=0.25)
plt.suptitle(f't = {time:.1f} hours')
plot_name = f'{plot_stem}_initial.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()

# ---------------------------------------------------------------------------- #
# FINAL PLOTTING
# ---------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 5))
time_idx = -1

for i, (field_name, colour_scheme, field_label, contour_to_remove, contours) in \
    enumerate(zip(
        final_field_names, final_colour_schemes,
        final_field_labels, final_contours_to_remove, final_contours)):

    # Make axes
    ax = fig.add_subplot(1, 1, 1+i, projection=projection)

    # Data extraction ----------------------------------------------------------
    field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
    coords_X, coords_Y = extract_gusto_coords(data_file, field_name)
    time = data_file['time'][time_idx] / (60.*60.)

    # Regrid onto lon-lat grid
    field_data = regrid_horizontal_slice(
        lon_2d, lat_2d, coords_X, coords_Y, field_data, periodic_fix='sphere'
    )

    # Plot data ----------------------------------------------------------------
    cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=contour_to_remove)
    cf, _ = plot_contoured_field(
        ax, lon_2d, lat_2d, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, projection=projection
    )

    add_colorbar_ax(ax, cf, field_label, location='bottom', cbar_labelpad=-10)
    tomplot_field_title(ax, None, minmax=True, field_data=field_data)

# Save figure ------------------------------------------------------------------
plt.suptitle(f't = {time:.1f} hours')
plot_name = f'{plot_stem}_final.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()
