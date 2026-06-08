"""
Plots the moist 3d bubble test case, plotting the theta_e field at a series of
times, in the x-z slice.
"""
from os.path import abspath, dirname
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_fig, tomplot_field_title, extract_gusto_vertical_slice
)

test = 'moist_3d_bubble'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_stem = f'{abspath(dirname(__file__))}/../../source/figures/{test}'
results_file_name = f'{abspath(dirname(__file__))}/../results/{test}/field_output.nc'

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
contours = np.linspace(319.5, 322.5, 13)
colour_scheme = 'OrRd'
field_label = r'$\theta_e$ (K)'
contour_method = 'contour'
xlims = [0, 10]
ylims = [0, 10]
field_name = 'Theta_e'
slice_along = 'y'
slice_at = 5.0625
time_idxs = [0, 1, 2]

# Things that are likely the same for all plots --------------------------------
set_tomplot_style(20)
data_file = Dataset(results_file_name, 'r')

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #
fig, axarray = plt.subplots(1, 3, figsize=(16, 5), sharex='all', sharey='all')

for i, (time_idx, ax) in enumerate(zip(time_idxs, axarray)):

    # Data extraction ----------------------------------------------------------
    field_data, coords_X, _, coords_Z = \
        extract_gusto_vertical_slice(
            data_file, field_name, time_idx=time_idx,
            slice_along=slice_along, slice_at=slice_at
        )
    time = data_file['time'][time_idx]

    # Plot data ----------------------------------------------------------------
    cmap, _ = tomplot_cmap(contours, colour_scheme)
    _, lines = tomplot_cmap(contours, colour_scheme, remove_contour=320.0)
    cf, lines = plot_contoured_field(
        ax, coords_X, coords_Z, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines
    )

    if i == 2:
        add_colorbar_fig(
            fig, cf, field_label, ax_idxs=[2], location='right',
            cbar_labelpad=-30, cbar_format='.1f'
        )

    tomplot_field_title(
        ax, r'$t=$'+f' ${time:.0f}$ s', minmax=False, field_data=field_data
    )

    # Labels -------------------------------------------------------------------
    ax.set_ylim(ylims)
    ax.set_yticks(ylims)

    if i == 0:
        ax.set_ylabel(r'$z$ (km)', labelpad=-20)
        ax.set_yticklabels(ylims)

    ax.set_xlabel(r'$x$ (km)', labelpad=-10)
    ax.set_xlim(xlims)
    ax.set_xticks(xlims)
    ax.set_xticklabels(xlims)

# Save figure ------------------------------------------------------------------
fig.suptitle(rf'Rising Thermal', y=0.99, fontsize=20)
fig.subplots_adjust(wspace=0.15)
plot_name = f'{plot_stem}.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight', dpi=300)
plt.close()
