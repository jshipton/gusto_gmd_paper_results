"""
Plots the baroclinic wave test case, after 10 days.

This plots:
(a) Surface Pressure , (b) Temperature
"""
from os.path import abspath, dirname
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
import numpy as np
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_ax, extract_gusto_coords, extract_gusto_field,
    reshape_gusto_data, area_restriction
)

test = 'dry_baroclinic_sphere'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
results_file_name = '/data/users/alex.brown/data/isambard/dry_baroclinic_sphere_dt450_C96_30L/field_output.nc'
plot_stem = f'{abspath(dirname(__file__))}/../figures/{test}'

# ---------------------------------------------------------------------------- #
# Final plot details
# ---------------------------------------------------------------------------- #
field_names = ['Pressure_Vt', 'Temperature']
slice_along = 'z'
colour_schemes = ['PiYG_r', 'RdBu_r']
centre_titles = ['Surface Pressure',  'Temperature']
field_labels = [r'$P$ (hPa)', r'$T$ (K)']

contours = [
    np.linspace(92000, 103000, 23),
    np.linspace(230, 330, 11)
]
line_contours = [
    np.linspace(92000, 103000, 12),
    np.linspace(230, 330, 11)
]
remove_lines = [100000.0, None]

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
contour_method = 'tricontour'
domain_limit = {'X' : (-180, 180), 'Y' : (-90, 90)}
time_idx = -1
# Things that are likely the same for all plots --------------------------------
set_tomplot_style(fontsize=24)
level = 0

y_ticks = [-90, 90]
y_tick_labels = ['90S', '90N']
xticks = [-180, 180]
xtick_labels = ['180W', '180E']

pressure_ticks = [92000, 94000, 96000, 98000, 100000, 102000]
pressure_labels = ['920', '940', '960', '980', '1000', '1020']
temp_ticks = [240, 280, 320]
temp_labels = ['240', '280', '320']

cbar_ticks = [pressure_ticks, temp_ticks]
cbar_ticklabels = [pressure_labels, temp_labels]

cmap_rescale_types = ['top', None]
cmap_rescalings = [0.72, None]

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #

data_file = Dataset(results_file_name, 'r')
time = data_file['time'][time_idx]
time_in_days = time / (24*60*60)

fig, axarray = plt.subplots(1, 2, figsize=(15, 6), sharex='col', sharey='row')

for i, (ax, field_name, colour_scheme, field_label, contour, line_contour,
        centre_title, cbar_tick, cbar_ticklabel, remove_contour, cmap_rescale_type,
        cmap_rescaling) in enumerate(zip(
        axarray, field_names, colour_schemes, field_labels, contours, line_contours,
        centre_titles, cbar_ticks, cbar_ticklabels, remove_lines,
        cmap_rescale_types, cmap_rescalings)):

    # Data extraction ------------------------------------------------------
    field_full = extract_gusto_field(data_file, field_name, time_idx)
    coords_X_full, coords_Y_full, coords_Z_full = \
        extract_gusto_coords(data_file, field_name)

    # Reshape
    field_full, coords_X_full, coords_Y_full, _ = \
        reshape_gusto_data(field_full, coords_X_full,
                            coords_Y_full, coords_Z_full)

    # Domain restriction
    field_data, coords_hori, coords_Z = \
        area_restriction(field_full[:, level], coords_X_full[:, level],
                                coords_Y_full[:, level], domain_limit)


    cmap, _ = tomplot_cmap(
        contour, colour_scheme, cmap_rescale_type=cmap_rescale_type,
        cmap_rescaling=cmap_rescaling
    )
    _, lines = tomplot_cmap(
        line_contour, colour_scheme, remove_contour=remove_contour
    )

    # Plot data ------------------------------------------------------------
    cf, _ = plot_contoured_field(
        ax, coords_hori, coords_Z, field_data, contour_method, contour,
        cmap=cmap, line_contours=lines
    )

    cb = add_colorbar_ax(ax, cf, field_label, location='bottom', cbar_ticks=cbar_tick, pad=0.07)
    cb.ax.set_xticklabels(cbar_ticklabel)

    # Labels ---------------------------------------------------------------
    # Adding titles
    ax.set_title(centre_title, fontsize='17.0')

    # Setting Y / X ticks
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize = '15.0')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, fontsize='15.0')

# Save figure ------------------------------------------------------------------
fig.suptitle(rf'Baroclinic Wave, $t=$ {int(time_in_days):02d} days', y=0.98, fontsize=24)
plot_name = f'{plot_stem}.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight', dpi=300)
plt.close()
