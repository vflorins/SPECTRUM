# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file = "main_test_parker_spiral_TrajectoryFieldline.lines"
plot_file = "plot_test_parker_spiral_TrajectoryFieldline.png"

# Import data
data = np.loadtxt(data_path + data_file, delimiter=",")

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.set_title("PARKER SPIRAL SOLAR WIND", fontsize=32)
ax.plot(data[:,0], data[:,1], color="tab:blue", linewidth=2)
ax.set_xlabel('$x$ ($au$)', fontsize=20)
ax.set_ylabel('$y$ ($au$)', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_aspect('equal')

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)