# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file = "main_test_dipole_drifts_TrajectoryGuiding.lines"
plot_file = "plot_test_dipole_drifts_TrajectoryGuiding.png"

# Import data
data = np.loadtxt(data_path + data_file, delimiter=",") / 4.26352e-5 # au -> RE
N = np.size(data,0)

# Find numbers
plot_label = ""
with open(r'benchmark_results.txt', 'r') as fp:
   lines = fp.readlines()
   for row in lines:
      word = "DIPOLE FIELD DRIFT PERIODS"  # String to search for
      if row.find(word) != -1:
         plot_label = lines[lines.index(row)+2]
         for l in range(lines.index(row)+4, lines.index(row)+9):
            plot_label = plot_label + lines[l]
plot_label = plot_label.strip()

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='3d')

ax.set_title("DIPOLE FIELD DRIFT PERIODS", fontsize=32)
ax.plot(data[:N//10,0], data[:N//10,1], data[:N//10,2], label=plot_label)
ax.set_xlabel('$x$ ($R_E$)', fontsize=20, labelpad=15)
ax.set_ylabel('$y$ ($R_E$)', fontsize=20, labelpad=15)
ax.set_zlabel('$z$ ($R_E$)', fontsize=20, labelpad=15)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='z', labelsize=20)
ax.set_aspect('equal')
ax.legend(loc = 1, fontsize=14)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)