import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import sys

args = sys.argv
numprocs = args[1]
print(args[1])
os.chdir("build")
cmd = "mpirun -np " + f"{numprocs}" + " ./task"
subprocess.run(cmd, shell=True, check=True)
df = pd.read_csv('data.csv', sep=' ', header=None)
df.drop([df.columns[100]], axis=1, inplace=True)
data = df

t_step = 1e-2
x_step = 1e-2

x = np.arange(0, 1, x_step)
t = np.arange(0, 1, t_step)

x, t = np.meshgrid(x, t)

fig = plt.figure(figsize = (12,8), dpi=512)
ax = plt.axes(projection='3d')
ax.set_xlabel('t', fontsize=11)
ax.set_ylabel('x', fontsize=11)
ax.set_zlabel('u(t, x)', fontsize=10)

ax.plot_surface(x, t, data, cmap = plt.cm.cividis)
plt.savefig('solution.png')

