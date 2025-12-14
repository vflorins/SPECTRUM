# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file1 = "main_test_dipole_visualization_Bx.dat"
data_file2 = "main_test_dipole_visualization_By.dat"
data_file3 = "main_test_dipole_visualization_Bz.dat"
plot_file = "plot_test_dipole_visualization_B.png"

# Import data
BX = np.transpose(np.loadtxt(data_path + data_file1))
BY = np.transpose(np.loadtxt(data_path + data_file2))
BZ = np.transpose(np.loadtxt(data_path + data_file3))
BM = np.sqrt(np.square(BX) + np.square(BY) + np.square(BZ))
Nx = np.size(BM, 1)
Nz = np.size(BM, 0)
X = np.linspace(-3.0, 3.0, num = Nx)
Z = np.linspace(-3.0, 3.0, num = Nz)
XX, ZZ = np.meshgrid(X, Z)

fig = plt.figure(figsize=(14, 12), layout='tight')
plt.suptitle("DIPOLE FIELD VISUALIZATION", fontsize=32)

ax1 = fig.add_subplot(221, projection='rectilinear')
ax1.set_xlabel("$x$ ($R_E$)", fontsize = 16)
ax1.set_ylabel("$z$ ($R_E$)", fontsize = 16)
ax1.set_title("$|B|$", fontsize = 16)
ax1.tick_params(labelsize=16)
hm = ax1.pcolormesh(XX, ZZ, BM, norm=LogNorm(vmin=1.0, vmax=100.0))
cb1 = fig.colorbar(hm, ax=ax1)
cb1.ax.tick_params(labelsize=16)

ax2 = fig.add_subplot(222, projection='rectilinear')
ax2.set_xlabel("$x$ ($R_E$)", fontsize = 16)
ax2.set_ylabel("$z$ ($R_E$)", fontsize = 16)
ax2.set_title("$B_x$", fontsize = 16)
ax2.tick_params(labelsize=16)
hm = ax2.pcolormesh(XX, ZZ, BX, vmin = -10.0, vmax = 10.0)
cb2 = fig.colorbar(hm, ax=ax2)
cb2.ax.tick_params(labelsize=16)

ax3 = fig.add_subplot(223, projection='rectilinear')
ax3.set_xlabel("$x$ ($R_E$)", fontsize = 16)
ax3.set_ylabel("$z$ ($R_E$)", fontsize = 16)
ax3.set_title("$B_y$", fontsize = 16)
ax3.tick_params(labelsize=16)
hm = ax3.pcolormesh(XX, ZZ, BY, vmin = -10.0, vmax = 10.0)
cb3 = fig.colorbar(hm, ax=ax3)
cb3.ax.tick_params(labelsize=16)

ax4 = fig.add_subplot(224, projection='rectilinear')
ax4.set_xlabel("$x$ ($R_E$)", fontsize = 16)
ax4.set_ylabel("$z$ ($R_E$)", fontsize = 16)
ax4.set_title("$B_z$", fontsize = 16)
ax4.tick_params(labelsize=16)
hm = ax4.pcolormesh(XX, ZZ, BZ, vmin = -10.0, vmax = 10.0)
cb4 = fig.colorbar(hm, ax=ax4)
cb4.ax.tick_params(labelsize=16)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)