"""
A tomplot example, for making a convergence plot.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
from os.path import abspath, dirname
from tomplot import set_tomplot_style, only_minmax_ticklabels
import pandas as pd


def format_dofs(n):
    """Format a DoF count as a LaTeX mantissa x 10^exponent string."""
    exponent = int(np.floor(np.log10(n)))
    mantissa = n / 10**exponent
    return rf'${mantissa:.1f} \times 10^{{{exponent}}}$'

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/../figures'

plot_name = f'{plot_dir}/gusto_scalability.png'
time_name = 'from_first_step_wall_time'

set_tomplot_style(20)
fig, axarray = plt.subplots(1, 2, figsize=(16, 7))

# ---------------------------------------------------------------------------- #
# Strong Scaling
# ---------------------------------------------------------------------------- #

csv_file = f'/data/users/thomas.bendall/gusto_gmd_paper_results/results/strong_scaling.csv'
xlabel = r"Number of Cores"
ylabel = r"Speed (Timesteps Per Wallclock Second)"

df = pd.read_csv(csv_file)
num_processes = df['nprocs'].values
num_steps = df['tmax'].values / df['dt'].values
speeds = num_steps / df[time_name].values

# Convert into a smaller list of unique process counts
unique_process_counts = sorted(set(num_processes))

cnum = df['ncells'].values
nlayers = df['nlayers'].values
total_cols = 6*cnum**2
total_edges = 12*cnum**2
u_dofs = total_edges*nlayers
w_dofs = total_cols*(nlayers+1)
p_dofs = total_cols*nlayers
total_dofs = u_dofs + 2*w_dofs + p_dofs

# Get mean, min and max speeds for each unique process count
mean_speeds = []
min_speeds = []
max_speeds = []
mean_total_dofs = []
for count in unique_process_counts:
    speeds_for_count = [speed for nproc, speed in zip(num_processes, speeds) if nproc == count]
    mean_speeds.append(sum(speeds_for_count) / len(speeds_for_count))
    min_speeds.append(min(speeds_for_count))
    max_speeds.append(max(speeds_for_count))
    dofs_for_count = [dofs for nproc, dofs in zip(num_processes, total_dofs) if nproc == count]
    mean_total_dofs.append(sum(dofs_for_count) / len(dofs_for_count))

# DoFs per processor for each unique process count (problem size is fixed)
dofs_per_proc = [dofs / count for dofs, count in zip(mean_total_dofs, unique_process_counts)]

# Plot mean speeds with error bars for min and max
axarray[0].errorbar(
    unique_process_counts, mean_speeds,
    yerr=[[mean - min_speed for mean, min_speed in zip(mean_speeds, min_speeds)],
          [max_speed - mean for mean, max_speed in zip(mean_speeds, max_speeds)]],
    fmt='+', color='black', ecolor='black', capsize=5
)

# Annotate each point with the DoFs per core, below the point
for i, (x, y, dofs) in enumerate(zip(unique_process_counts, mean_speeds, dofs_per_proc)):
    if i == 0:
        ha = 'left'
        va = 'top'
        offset = (0, -18)
    elif i == len(unique_process_counts) - 1:
        ha = 'right'
        va = 'bottom'
        offset = (8, 10)
    else:
        ha = 'center'
        va = 'top'
        offset = (0, -18)
    axarray[0].annotate(
        f'DoFs per core:\n{format_dofs(dofs)}', (x, y), textcoords='offset points',
        xytext=offset, ha=ha, va=va, fontsize=11
    )

# Ideal scaling line, going through first point
ideal_speeds = [mean_speeds[0] * n / unique_process_counts[0] for n in unique_process_counts]
axarray[0].plot(unique_process_counts, ideal_speeds, linestyle='--', color='gray')

axarray[0].set_xscale('log', base=2)
axarray[0].set_yscale('log')
axarray[0].xaxis.set_major_formatter(ScalarFormatter())
axarray[0].xaxis.set_minor_formatter(ScalarFormatter())
# Pad title
axarray[0].set_ylim(1e-2, 1e0)

axarray[0].set_title('Strong Scalability (fixed problem size)', pad=10)
axarray[0].set_xlabel(xlabel)
axarray[0].set_ylabel(ylabel, labelpad=5)
axarray[0].text(
    0.5, 0.95, f'Total DoFs: {format_dofs(mean_total_dofs[0])}',
    transform=axarray[0].transAxes, ha='center', va='top', fontsize=16
)
axarray[0].grid()

# ---------------------------------------------------------------------------- #
# Weak Scaling
# ---------------------------------------------------------------------------- #

csv_file = f'/data/users/thomas.bendall/gusto_gmd_paper_results/results/weak_scaling.csv'
xlabel = r"Number of Cores"
ylabel = r"Speed (Timesteps Per Wallclock Second)"

df = pd.read_csv(csv_file)
num_processes = df['nprocs'].values
num_steps = df['tmax'].values / df['dt'].values
speeds = num_steps / df[time_name].values

cnum = df['ncells'].values
nlayers = df['nlayers'].values
total_cols = 6*cnum**2
total_edges = 12*cnum**2
u_dofs = total_edges*nlayers
w_dofs = total_cols*(nlayers+1)
p_dofs = total_cols*nlayers
total_dofs = u_dofs + 2*w_dofs + p_dofs

# Convert into a smaller list of unique process counts
unique_process_counts = sorted(set(num_processes))

# Get mean, min and max speeds for each unique process count
mean_speeds = []
min_speeds = []
max_speeds = []
mean_total_dofs = []
for count in unique_process_counts:
    speeds_for_count = [speed for nproc, speed in zip(num_processes, speeds) if nproc == count]
    mean_speeds.append(sum(speeds_for_count) / len(speeds_for_count))
    min_speeds.append(min(speeds_for_count))
    max_speeds.append(max(speeds_for_count))
    dofs_for_count = [dofs for nproc, dofs in zip(num_processes, total_dofs) if nproc == count]
    mean_total_dofs.append(sum(dofs_for_count) / len(dofs_for_count))

# DoFs per processor (fixed local volume per process, so this should be constant)
dofs_per_proc = [dofs / count for dofs, count in zip(mean_total_dofs, unique_process_counts)]

# Plot mean speeds with error bars for min and max
axarray[1].errorbar(
    unique_process_counts, mean_speeds,
    yerr=[[mean - min_speed for mean, min_speed in zip(mean_speeds, min_speeds)],
          [max_speed - mean for mean, max_speed in zip(mean_speeds, max_speeds)]],
    fmt='+', color='black', ecolor='black', capsize=5
)

# Annotate each point with the total number of DoFs, below the point
for i, (x, y, dofs) in enumerate(zip(unique_process_counts, mean_speeds, mean_total_dofs)):
    if i == 0:
        ha = 'left'
    elif i == len(unique_process_counts) - 1:
        ha = 'right'
    else:
        ha = 'center'
    axarray[1].annotate(
        f'Total DoFs:\n{format_dofs(dofs)}', (x, y), textcoords='offset points',
        xytext=(0, -18), ha=ha, va='top', fontsize=11
    )

axarray[1].set_xscale('log', base=2)
axarray[1].xaxis.set_major_formatter(ScalarFormatter())
axarray[1].xaxis.set_minor_formatter(ScalarFormatter())
axarray[1].set_title('Weak Scalability (fixed local volume per process)', pad=10)
axarray[1].set_xlabel(xlabel)
axarray[1].set_ylabel(ylabel, labelpad=-20)
axarray[1].text(
    0.5, 0.95, f'DoFs per core: {format_dofs(sum(dofs_per_proc) / len(dofs_per_proc))}',
    transform=axarray[1].transAxes, ha='center', va='top', fontsize=16
)
axarray[1].set_ylim(0, 0.3)
axarray[1].set_yticks([0, 0.1, 0.2, 0.35])
axarray[1].set_yticklabels([r'$0.0$', '', '', r'$0.3$'])
axarray[1].grid()

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
fig.subplots_adjust(wspace=0.3)  # Adjust horizontal space between subplots
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


