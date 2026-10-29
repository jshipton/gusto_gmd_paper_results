"""
Plots the Straka bubble test case, comparing order 0 and order 1 elements.

This plots the theta perturbation @ t = 900 s:
(a) 0-th order elements, (b) 1-st order elements
"""
from os.path import abspath, dirname
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from tomplot import (
    set_tomplot_style, tomplot_cmap, plot_contoured_field,
    add_colorbar_fig, tomplot_field_title, extract_gusto_coords,
    extract_gusto_field
)

test = 'straka'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_stem = f'{abspath(dirname(__file__))}/../figures/{test}'

# ---------------------------------------------------------------------------- #
# Plot details
# ---------------------------------------------------------------------------- #
results_opts = ['degree_0', 'degree_1']
titles = ['Degree 0 Elements', 'Degree 1 Elements']
cbars = [False, True]

# ---------------------------------------------------------------------------- #
# General options
# ---------------------------------------------------------------------------- #
contours = np.linspace(-7.5, 0.5, 17)
colour_scheme = 'Blues_r'
field_label = r'$\Delta \theta$ (K)'
contour_method = 'tricontour'
xlims = [0, 12]
ylims = [0, 5]
field_name = 'theta_perturbation'
time_idx = -1

# Things that are likely the same for all plots --------------------------------
set_tomplot_style(24)

# ---------------------------------------------------------------------------- #
# PLOTTING
# ---------------------------------------------------------------------------- #
fig, axarray = plt.subplots(2, 1, figsize=(12, 10), sharex='all', sharey='all')

for i, (ax, results_opt, cbar, title) in \
        enumerate(zip(axarray.flatten(), results_opts, cbars, titles)):

    results_file_name = f'{abspath(dirname(__file__))}/../results/{test}_{results_opt}/field_output.nc'
    data_file = Dataset(results_file_name, 'r')

    # Data extraction ----------------------------------------------------------
    field_data = extract_gusto_field(data_file, field_name, time_idx=time_idx)
    coords_X, coords_Y = extract_gusto_coords(data_file, field_name)
    time = data_file['time'][time_idx]

    # Translate coordinates
    coords_X -= 25.6

    # Plot data ----------------------------------------------------------------
    cmap, _ = tomplot_cmap(contours, colour_scheme)
    _, lines = tomplot_cmap(contours, colour_scheme, remove_contour=0.0)
    cf, lines = plot_contoured_field(
        ax, coords_X, coords_Y, field_data, contour_method, contours,
        cmap=cmap, line_contours=lines, negative_linestyles='solid'
    )

    if cbar:
        add_colorbar_fig(
            fig, cf, field_label, ax_idxs=[1], location='bottom',
            cbar_labelpad=-10, cbar_padding=0.02, cbar_format='.1f'
        )
    tomplot_field_title(
        ax, title, minmax=False, field_data=field_data
    )

    # Labels -------------------------------------------------------------------
    ax.set_ylabel(r'$z$ (km)', labelpad=-10)
    ax.set_ylim(ylims)
    ax.set_yticks(ylims)
    ax.set_yticklabels(ylims)

    if i == 1:
        ax.set_xlabel(r'$x$ (km)', labelpad=-10)
        ax.set_xlim(xlims)
        ax.set_xticks(xlims)
        ax.set_xticklabels(xlims)

# Save figure ------------------------------------------------------------------
fig.suptitle(rf'Density Current, $t=$ {time:.0f} s', y=0.97, fontsize=28)
fig.subplots_adjust(hspace=0.15)
plot_name = f'{plot_stem}.png'
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight', dpi=300)
plt.close()
