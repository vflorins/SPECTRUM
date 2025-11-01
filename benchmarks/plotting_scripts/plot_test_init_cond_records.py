# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file1 = "main_test_init_cond_records_TrajectoryFocused_init_time.dat"
data_file2 = "main_test_init_cond_records_TrajectoryFocused_pos_records.dat"
plot_file = "plot_test_init_cond_records_TrajectoryFocused.png"

# Import data
data1 = np.loadtxt(data_path + data_file1, skiprows=2)
data2 = np.loadtxt(data_path + data_file2, skiprows=2)

# Process data
n_bins = 100
bin_size = 2.0 / n_bins
S = np.linspace(-1.0, 1.0, num=n_bins+1)
s = 0.5 * (S[1:] + S[:-1])
Nx = s[:] * 0.0
Ny = s[:] * 0.0
Nz = s[:] * 0.0

for i in range(np.size(data2, 0)):
   idx = np.floor((data2[i,0] + 1.0) / bin_size).astype(int)
   Nx[idx] = Nx[idx] + 1
   idx = np.floor((data2[i,1] + 1.0) / bin_size).astype(int)
   Ny[idx] = Ny[idx] + 1
   idx = np.floor((data2[i,2] + 1.0) / bin_size).astype(int)
   Nz[idx] = Nz[idx] + 1

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')
plt.suptitle("INITIAL CONDITION RECORDS", fontsize=32)

ax1 = fig.add_subplot(221, projection='rectilinear')

ax1.plot(data1[:,0], data1[:,2], color="tab:red", linestyle='-')
ax1.set_xlabel('$t$', fontsize=20)
ax1.set_ylabel('$N(t)$', fontsize=20)
ax1.set_ylim(bottom=0.0)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)

ax2 = fig.add_subplot(222, projection='rectilinear')

ax2.plot(s, Nx, color="tab:red", linestyle='-')
ax2.set_xlabel('$x$', fontsize=20)
ax2.set_ylabel('$N(x)$', fontsize=20)
ax2.set_ylim(bottom=0.0)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)

ax3 = fig.add_subplot(223, projection='rectilinear')

ax3.plot(s, Ny, color="tab:red", linestyle='-')
ax3.set_xlabel('$y$', fontsize=20)
ax3.set_ylabel('$N(y)$', fontsize=20)
ax3.set_ylim(bottom=0.0)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)

ax4 = fig.add_subplot(224, projection='rectilinear')

ax4.plot(s, Nz, color="tab:red", linestyle='-')
ax4.set_xlabel('$z$', fontsize=20)
ax4.set_ylabel('$N(z)$', fontsize=20)
ax4.set_ylim(bottom=0.0)
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)