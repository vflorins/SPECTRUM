# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file = "main_test_full_diff_TrajectoryParker_cumulative_distro2.dat"
plot_file = "plot_test_full_diff_TrajectoryParker_cumulative_distro2.png"

# Import data
data = np.loadtxt(data_path + data_file)

# Find numbers
line_with_number = ""
slope = 0
with open(r'benchmark_results.txt', 'r') as fp:
   lines = fp.readlines()
   for row in lines:
      word = "FULL DIFFUSION"  # String to search for
      if row.find(word) != -1:
         line_with_number = lines[lines.index(row)+3].split()
         slope = float(line_with_number[2])
theory = 2.0 * slope * data[:,0]

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.set_title("FULL DIFFUSION", fontsize=32)
ax.plot(data[:,0], theory, color="black", linestyle='-',
        linewidth = 4, label="theoretical")
ax.plot(data[:,0], data[:,1], color="tab:red", linestyle='--',
        linewidth = 4, label="simulation x")
ax.plot(data[:,0], data[:,5], color="tab:blue", linestyle='-.',
        linewidth = 4, label="simulation y")
ax.plot(data[:,0], data[:,9], color="tab:green", linestyle=':',
        linewidth = 4, label="simulation z")
ax.set_xlabel('$t$', fontsize=20)
ax.set_ylabel('$\\langle r^2 \\rangle$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=20)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)