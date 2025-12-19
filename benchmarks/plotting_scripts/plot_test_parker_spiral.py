# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file1 = "main_test_solarwind_parker_spiral_TrajectoryFieldline.lines"
data_file2 = "main_test_cartesian_parker_spiral_TrajectoryFieldline.lines"
plot_file = "plot_test_parker_spiral_TrajectoryFieldline.png"

# Import data
data1 = np.loadtxt(data_path + data_file1, delimiter=",") / 1.496e+13 # au
data2 = np.loadtxt(data_path + data_file2, delimiter=",") / 1.496e+13 # au

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.set_title("PARKER SPIRAL SOLAR WIND", fontsize=32)
ax.plot(data1[:,0], data1[:,1], color="tab:blue",
        linestyle="-", linewidth=2, label="solarwind")
ax.plot(data2[:,0], data2[:,1], color="tab:red",
        linestyle="--", linewidth=2, label="cartesian")
ax.set_xlabel('$x$ ($au$)', fontsize=20)
ax.set_ylabel('$y$ ($au$)', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_aspect('equal')
ax.legend(fontsize=16)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)