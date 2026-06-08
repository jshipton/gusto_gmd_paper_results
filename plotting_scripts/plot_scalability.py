"""
A tomplot example, for making a convergence plot.
"""

import matplotlib.pyplot as plt
from os.path import abspath, dirname
from tomplot import set_tomplot_style, only_minmax_ticklabels
import pandas as pd

# ---------------------------------------------------------------------------- #
# Directory for results and plots
# ---------------------------------------------------------------------------- #
plot_dir = f'{abspath(dirname(__file__))}/../../source/figures'
plot_name = f'{plot_dir}/gusto_scalability.png'

set_tomplot_style(20)
fig, axarray = plt.subplots(1, 2, figsize=(16, 7))

# ---------------------------------------------------------------------------- #
# Strong Scaling
# ---------------------------------------------------------------------------- #

csv_file = '/data/users/alex.brown/data/isambard/baroclinic_sphere_scaling/siqn_dryrun_strong_scaling.csv'
xlabel = r"Number of processes"
ylabel = r"Speed (Simulated Day Per Second)"

df = pd.read_csv(csv_file)
num_processes = df['nprocs'].values
speeds = 1.0 / df['wall_time'].values

# Convert into a smaller list of unique process counts
unique_process_counts = sorted(set(num_processes))

# Get mean, min and max speeds for each unique process count
mean_speeds = []
min_speeds = []
max_speeds = []
for count in unique_process_counts:
    speeds_for_count = [speed for nproc, speed in zip(num_processes, speeds) if nproc == count]
    mean_speeds.append(sum(speeds_for_count) / len(speeds_for_count))
    min_speeds.append(min(speeds_for_count))
    max_speeds.append(max(speeds_for_count))

# Plot mean speeds with error bars for min and max
axarray[0].errorbar(
    unique_process_counts, mean_speeds,
    yerr=[[mean - min_speed for mean, min_speed in zip(mean_speeds, min_speeds)],
          [max_speed - mean for mean, max_speed in zip(mean_speeds, max_speeds)]],
    fmt='+', color='black', ecolor='black', capsize=5
)

# Ideal scaling line, going through first point
ideal_speeds = [mean_speeds[0] * n / unique_process_counts[0] for n in unique_process_counts]
axarray[0].plot(unique_process_counts, ideal_speeds, linestyle='--', color='gray')

axarray[0].set_xscale('log', base=2)
axarray[0].set_yscale('log')
# Pad title
axarray[0].set_title('Strong Scalability (fixed problem size)', pad=10)
axarray[0].set_xlabel(xlabel)
axarray[0].set_ylabel(ylabel, labelpad=5)
axarray[0].grid()

# ---------------------------------------------------------------------------- #
# Weak Scaling
# ---------------------------------------------------------------------------- #

csv_file = '/data/users/alex.brown/data/isambard/baroclinic_sphere_scaling/siqn_weak_scaling.csv'
xlabel = r"Number of processes"
ylabel = r"Speed (Simulated Day Per Second)"

df = pd.read_csv(csv_file)
num_processes = df['nprocs'].values
speeds = 1.0 / df['wall_time'].values

# Convert into a smaller list of unique process counts
unique_process_counts = sorted(set(num_processes))

# Get mean, min and max speeds for each unique process count
mean_speeds = []
min_speeds = []
max_speeds = []
for count in unique_process_counts:
    speeds_for_count = [speed for nproc, speed in zip(num_processes, speeds) if nproc == count]
    mean_speeds.append(sum(speeds_for_count) / len(speeds_for_count))
    min_speeds.append(min(speeds_for_count))
    max_speeds.append(max(speeds_for_count))

# Plot mean speeds with error bars for min and max
axarray[1].errorbar(
    unique_process_counts, mean_speeds,
    yerr=[[mean - min_speed for mean, min_speed in zip(mean_speeds, min_speeds)],
          [max_speed - mean for mean, max_speed in zip(mean_speeds, max_speeds)]],
    fmt='+', color='black', ecolor='black', capsize=5
)

axarray[1].set_xscale('log', base=2)
axarray[1].set_title('Weak Scalability (fixed local volume per process)', pad=10)
axarray[1].set_xlabel(xlabel)
axarray[1].set_ylim(2e-4, 3.5e-4)
axarray[1].set_yticks([2e-4, 2.5e-4, 3e-4, 3.5e-4])
axarray[1].set_yticklabels([r'$2\times 10^{-4}$', None, None, r'$3.5\times 10^{-4}$'])
axarray[1].grid()

# ---------------------------------------------------------------------------- #
# Save figure
# ---------------------------------------------------------------------------- #
fig.subplots_adjust(wspace=0.3)  # Adjust horizontal space between subplots
print(f'Saving figure to {plot_name}')
fig.savefig(plot_name, bbox_inches='tight')
plt.close()


