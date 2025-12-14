# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

# Grab command line arguments
data_path = sys.argv[1]
plot_path = sys.argv[2]
data_file1 = "main_test_pa_distro_isotrop_TrajectoryGuidingScatt_initial_distro.dat"
data_file2 = "main_test_pa_distro_isotrop_TrajectoryGuidingScatt_timemark1_distro.dat"
data_file3 = "main_test_pa_distro_isotrop_TrajectoryGuidingScatt_timemark2_distro.dat"
data_file4 = "main_test_pa_distro_isotrop_TrajectoryGuidingScatt_timemark3_distro.dat"
data_file5 = "main_test_pa_distro_isotrop_TrajectoryGuidingScatt_final_distro.dat"
plot_file = "plot_test_pa_distro_isotrop_TrajectoryGuidingScatt.png"

# Import data
data1 = np.loadtxt(data_path + data_file1, skiprows=2)
data2 = np.loadtxt(data_path + data_file2, skiprows=2)
data3 = np.loadtxt(data_path + data_file3, skiprows=2)
data4 = np.loadtxt(data_path + data_file4, skiprows=2)
data5 = np.loadtxt(data_path + data_file5, skiprows=2)
N = 100000 * 2 / np.size(data1, 0)

# Find numbers
line_with_number = ""
scatt_time = 1.0
time1 = 0.0
time2 = 0.0
time3 = 0.0
time4 = 0.0
time5 = 0.0
with open(r'benchmark_results.txt', 'r') as fp:
   lines = fp.readlines()
   for row in lines:
      word = "PITCH ANGLE DISTRIBUTION ISOTROPIZATION"  # String to search for
      if row.find(word) != -1:
         line_with_number = lines[lines.index(row)+3].split()
         scatt_time = 1.0 / float(line_with_number[2])
         line_with_number = lines[lines.index(row)+5].split()
         time1 = float(line_with_number[2])
         line_with_number = lines[lines.index(row)+6].split()
         time2 = float(line_with_number[2])
         line_with_number = lines[lines.index(row)+7].split()
         time3 = float(line_with_number[2])
         line_with_number = lines[lines.index(row)+8].split()
         time4 = float(line_with_number[2])
         line_with_number = lines[lines.index(row)+9].split()
         time5 = float(line_with_number[2])

# Make figure
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.set_title("PITCH ANGLE DISTRIBUTION ISOTROPIZATION", fontsize=32)
ax.plot(data1[:,0], data1[:,2] / N, color="black", linestyle='-',
        linewidth = 3, label="t = {:.2e} s".format(time1))
ax.plot(data2[:,0], data2[:,2] / N, color="tab:red", linestyle='-',
        linewidth = 3, label="t = {:.2e} s".format(time2))
ax.plot(data3[:,0], data3[:,2] / N, color="tab:blue", linestyle='-',
        linewidth = 3, label="t = {:.2e} s".format(time3))
ax.plot(data4[:,0], data4[:,2] / N, color="tab:green", linestyle='-',
        linewidth = 3, label="t = {:.2e} s".format(time4))
ax.plot(data5[:,0], data5[:,2] / N, color="tab:orange", linestyle='-',
        linewidth = 3, label="t = {:.2e} s".format(time5))
ax.set_xlabel('$\\mu$', fontsize=20)
ax.set_ylabel('$f(\\mu)$ (normalized)', fontsize=20)
ax.set_xlim(-1.0,1.0)
ax.set_ylim(0.0,1.0)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=20)

# Save figure
plt.savefig(plot_path + plot_file)
plt.close(fig)