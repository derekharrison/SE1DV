import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.animation as manimation
import time
import sys

start_time = time.time()

plt.rcParams['animation.ffmpeg_path']='C:/Users/d-w-h/Downloads/ffmpeg-20200818-1c7e55d-win64-static/ffmpeg-20200818-1c7e55d-win64-static/bin/ffmpeg.exe'
writer=manimation.FFMpegWriter(bitrate=20000, fps=15)

fig = plt.figure(figsize=(8,8))

file_timesteps = 'number_of_timesteps.txt'
Nt = np.genfromtxt(file_timesteps, unpack=True)

file_limits = 'limits.txt'
max_real, min_real, max_im, min_im, L = np.genfromtxt(file_limits, unpack=True)

max_val = 0
min_val = 0
if(max_real > max_im):
    max_val = max_real
else:
    max_val = max_im

if(min_real < min_im):
    min_val = min_real
else:
    min_val = min_im

def animate(i):
    my_file = 'psi_vs_t_' + str(i) + '.txt'
    print(i)
    fig.clear()
    x_c, psi_real, psi_im = np.genfromtxt(my_file, unpack=True)
    ax = plt.axes(xlim=(-0.5*L, 0.5*L), ylim=(min_im, max_im))
    cont = plt.plot(x_c, psi_real)
    cont = plt.plot(x_c, psi_im)
    
    return cont

size_t = int(Nt)
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('SE_solution.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
