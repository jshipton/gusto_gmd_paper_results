"""
Plots Gusto spatial and temporal convergence results.
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from os.path import abspath, dirname
import numpy as np
from tomplot import set_tomplot_style, plot_convergence

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/../figures'
plot_name = f'{plot_dir}/gusto_convergence.png'
base_results_dir = f'{abspath(dirname(__file__))}/../results/'
base_results_name = 'linear_sw_gw'

set_tomplot_style()
fig, axarray = plt.subplots(1, 2, figsize=(12, 5))

# ---------------------------------------------------------------------------- #
# Spatial convergence
# ---------------------------------------------------------------------------- #
log_by = 'data'
xlabel = r"$\log(\Delta x)$"
ylabel = r"$\log(||D - D_0||)$"
legend_loc = 'upper center'
xlims = [12.7, 13.5]  # Set manually
ylims = [11, 22]

elements = ['RTCF0', 'RTCF1', 'RTCF2', 'BDM1', 'BDM2']
colours = ['red', 'blue', 'purple', 'orange', 'green', 'brown', 'pink']
markers = ['s', 'o', '^', 'D', 'v', 'P', '*']
labels = [r'RTCF$_0$:', r'RTCF$_1$:', r'RTCF$_2$:', r'BDM$_1$:', r'BDM$_2$:']

for element, colour, marker, label in zip(elements, colours, markers, labels):
    error_data = []

    if element in ['BDM1', 'BDM2']:
        dt = 60
        n_values = [11, 12, 13, 14, 15, 16, 17]
        dx_values = [2*np.pi*6371220.*np.cos(np.atan2(1,2))/(5*n) for n in n_values]

    elif element in ['RTCF0', 'RTCF1', 'RTCF2']:
        dt = 60
        n_values = [16, 17, 18, 19, 20, 21, 22, 23]
        dx_values = [2*np.pi*6371220./(4*n) for n in n_values]

    for n_value in n_values:
        results_name = f'{base_results_dir}/{base_results_name}_trapezium_rule_{element}_n{n_value}_dt{dt}/error_diagnostic.nc'
        results_file = Dataset(results_name, 'r')
        data_point = results_file['D_error'][-1]
        error_data.append(data_point)

    plot_convergence(
        axarray[0], dx_values, error_data, label=label,
        color=colour, marker=marker, log_by=log_by
    )

axarray[0].set_title('Spatial convergence')
axarray[0].set_xlabel(xlabel, labelpad=-10)
axarray[0].set_ylabel(ylabel, labelpad=-10)
axarray[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize='small')

# Manually set ticks
axarray[0].set_xticks([13.0, 13.1, 13.2, 13.3, 13.4])
axarray[0].set_xticklabels(['13.0', None, None, None, '13.4'])
axarray[0].set_yticks([12, 14, 16, 18, 20])
axarray[0].set_yticklabels(['12', None, None, None, '20'])
axarray[0].grid()

# ---------------------------------------------------------------------------- #
# Temporal convergence
# ---------------------------------------------------------------------------- #

log_by = 'data'
xlabel = r"$\log(\Delta t)$"
ylabel = r"$\log(||u - u_0||)$"
legend_loc = 'upper center'
base_results_name = 'w6'

schemes = ['backward_euler', 'siqn', 'heun', 'ssprk3', 'imex_ssp3', 'imex_ark4']
colours = ['red', 'blue', 'purple', 'orange', 'green', 'brown', 'pink']
markers = ['s', 'o', '^', 'D', 'v', 'P', '*']
labels = [r'Backward Euler:', r'SIQN:', r'Heun:', r'SSPRK3:', r'IMEX SSP3:', r'IMEX ARK4:']
explicit_dt_values = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]
implicit_dt_values = [240.0, 120.0, 60.0, 30.0]
all_dt_values = [
    implicit_dt_values,
    implicit_dt_values,
    explicit_dt_values,
    explicit_dt_values,
    implicit_dt_values,
    implicit_dt_values
]

dt_str_dict = {
    240.0: '240',
    120.0: '120',
    60.0: '60',
    30.0: '30',
    15.0: '15',
    7.5: '7p5',
    3.75: '3p75',
    1.875: '1p875'
}

for i, (scheme, dt_values, colour, marker, label) in enumerate(zip(schemes, all_dt_values, colours, markers, labels)):

    error_data = []
    actual_dt_values = []

    for dt_value in dt_values:
        try:
            dt_str = dt_str_dict[dt_value]
            results_name = f'{base_results_dir}/{base_results_name}_{scheme}_RTCF1_n128_dt{dt_str}/error_diagnostic.nc'
            results_file = Dataset(results_name, 'r')
            data_point = results_file['u_error'][:][-1]
            error_data.append(data_point)
            actual_dt_values.append(dt_value)
        except FileNotFoundError:
            print(f"Warning: File not found for {scheme} with dt={dt_value}. Skipping this point.")

    plot_convergence(
        axarray[1], actual_dt_values, error_data, label=label,
        color=colour, marker=marker, log_by=log_by
    )

axarray[1].set_title('Temporal convergence')
axarray[1].set_xlabel(xlabel, labelpad=-10)
axarray[1].set_ylabel(ylabel, labelpad=-10)
axarray[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize='small')

# Manually set ticks
axarray[1].set_xticks([1, 2, 3, 4, 5])
axarray[1].set_xticklabels(['1', None, None, None, 5])
axarray[1].set_yticks([0, 5, 10, 15])
axarray[1].set_yticklabels(['0', None, None, '15'])
axarray[1].grid()

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


