"""
Plot from the Terminator Toy test case.
"""
from os.path import abspath, dirname
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field, add_colorbar_fig,
    tomplot_field_title, extract_gusto_coords, extract_gusto_field,
)

test = 'terminator_toy'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
results_file_name = f'{abspath(dirname(__file__))}/../results/{test}/field_output.nc'
plot_stem = f'{abspath(dirname(__file__))}/../../source/figures/{test}'

# ---------------------------------------------------------------------------- #
# Plot details
# ---------------------------------------------------------------------------- #
subplots_x = 2
subplots_y = 3
field_names = ['Y', 'Y',
               'Y2', 'Y2',
               'Y_plus_Y2_plus_Y2', 'Y_plus_Y2_plus_Y2']
time_idxs = [0, -1,
             0, -1,
             0, -1]
cbars = [False, True,
         False, True,
         False, True]

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #

projection = ccrs.Robinson()

Xt_contours = np.linspace(-1e-13, 1e-13, 11)
Xt_colour_scheme = 'BrBG'
Xt_field_label = r'$\Delta m_{X_T}$ (kg kg$^{-1}$)'

X_contours = np.linspace(0.0, 4.e-6, 17)
X_colour_scheme = 'PuRd'
X_field_label = r'$m_X$ (kg kg$^{-1}$)'

X2_contours_first = np.linspace(0.0, 2.000000000001e-6, 17)
X2_contours = np.linspace(0.0, 2.0e-6, 17)
X2_colour_scheme = 'PuBu'
X2_field_label = r'$m_{X_2}$ (kg kg$^{-1}$)'

contour_method = 'tricontour'

# Things that are likely the same for all plots --------------------------------
set_tomplot_style()
data_file = Dataset(results_file_name, 'r')

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #

fig = plt.figure(figsize=(12, 12))

# Save filled contour data for adding colourbars at end
all_cf = []

for i, (time_idx, field_name) in enumerate(zip(time_idxs, field_names)):

    ax = fig.add_subplot(subplots_y, subplots_x, 1+i, projection=projection)

    if time_idx == 'midpoint':
        time_idx = int((len(data_file['time'][:]) - 1) / 2)

    # Data extraction ----------------------------------------------------------
    field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
    coords_X, coords_Y = extract_gusto_coords(data_file, field_name)
    # Quote time in days:
    time = data_file['time'][time_idx] / (24.*60.*60.)

    # Select options for each field --------------------------------------------
    if field_name == 'Y_plus_Y2_plus_Y2':
        field_data -= 4.0e-6
        contours = Xt_contours
        colour_scheme = Xt_colour_scheme
        cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=None)
        data_format = '.2e'

    elif field_name == 'Y':
        contours = X_contours
        colour_scheme = X_colour_scheme
        cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=None)
        data_format = '.2e'

    elif field_name == 'Y2':
        if time_idx == 0:
            contours = X2_contours_first
        else:
            contours = X2_contours
        colour_scheme = X2_colour_scheme
        cmap, lines = tomplot_cmap(contours, colour_scheme, remove_contour=None)
        data_format = '.2g'

    # Plot data ----------------------------------------------------------------

    cf, _ = plot_contoured_field(
        ax, coords_X, coords_Y, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, projection=projection
    )
    all_cf.append(cf)

    tomplot_field_title(
        ax, None, minmax=True, minmax_format=data_format,
        field_data=field_data
    )

    if i in [0, 1]:
        ax.text(
            0.5, 1.2, r'$t =$ '+f'{time:.0f} days',
            horizontalalignment='center', transform=ax.transAxes
        )

# Colour bars need adding at end when the full figure exists -------------------
fig.subplots_adjust(wspace=0.12, hspace=-0.05)

for i, (cbar, field_name) in enumerate(zip(cbars, field_names)):

    cbar_labelpad = -50
    data_format = '.1e'

    # Get information for field
    if field_name == 'Y_plus_Y2_plus_Y2':
        field_label = Xt_field_label

    elif field_name == 'Y':
        field_label = X_field_label

    elif field_name == 'Y2':
        field_label = X2_field_label

    if cbar:
        add_colorbar_fig(
            fig, all_cf[i], field_label, ax_idxs=[i-1, i], location='right',
            cbar_labelpad=cbar_labelpad, cbar_format=data_format
        )

# Save figure ------------------------------------------------------------------
fig.suptitle(rf'Terminator Toy', y=0.92, fontsize=20)
plot_name = f'{plot_stem}.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight', dpi=300)
plt.close()
