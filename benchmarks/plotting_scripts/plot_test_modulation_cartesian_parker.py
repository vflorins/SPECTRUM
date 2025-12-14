# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file = "main_test_modulation_cartesian_TrajectoryParker_spectrum_pp.dat"
plot_file = "plot_test_modulation_cartesian_TrajectoryParker_spectrum_pp.png"

# Import data
data = np.loadtxt(data_path + data_file)

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.set_title("MODULATION CARTESIAN PARKER FIELD", fontsize=32)
ax.loglog(data[:,0], data[:,1], color="black", linestyle='-',
          linewidth = 3, marker = '', label="unmodulated")
ax.loglog(data[:,0], data[:,2], color="tab:red", linestyle='--',
          linewidth = 3, marker = '', label="modulated theory")
ax.loglog(data[:,0], data[:,3], color="tab:blue", linestyle='',
          linewidth = 3, marker = 'o', label="modulated simulation")
ax.set_xlabel('$E$ (MeV)', fontsize=20)
ax.set_ylabel('$4\\pi fp^2$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=20)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)