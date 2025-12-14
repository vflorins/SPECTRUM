# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file = "main_test_turb_waves_TrajectoryFieldline.lines"
plot_file = "plot_test_turb_waves_TrajectoryFieldline.png"

# Import data
data = np.loadtxt(data_path + data_file, delimiter=",")
N = np.size(data, 0)

# Find numbers
plot_label = ""
with open(r'benchmark_results.txt', 'r') as fp:
   lines = fp.readlines()
   for row in lines:
      word = "TURBULENCE VIA SUPERPOSITION OF WAVES"  # String to search for
      if row.find(word) != -1:
         for l in range(lines.index(row)+2, lines.index(row)+8):
            plot_label = plot_label + lines[l]
plot_label = plot_label.strip()

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='3d')

ax.set_title("TURBULENCE VIA SUPERPOSITION OF WAVES", fontsize=32)
ax.plot(data[:N,0], data[:N,1], data[:N,2], label=plot_label)
ax.set_xlabel('$x$', fontsize=20)
ax.set_ylabel('$y$', fontsize=20)
ax.set_zlabel('$z$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='z', labelsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_aspect('auto')
ax.legend(loc = 1, fontsize=14)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)