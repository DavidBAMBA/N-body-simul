import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  
import pandas as pd

num_steps = 1200  

def load_particle_positions(step):
    filename = f'positions_{step}.csv'
    data = pd.read_csv(filename)
    x_positions = data['X'].values
    y_positions = data['Y'].values
    z_positions = data['Z'].values
    return x_positions, y_positions, z_positions

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])
ax.set_zlim([-100, 100])
particles, = ax.plot([], [], [], 'bo', ms=1)  
def init():
    particles.set_data([], [])
    particles.set_3d_properties([])  
    return particles,

def animate(step):
    x, y, z = load_particle_positions(step)
    particles.set_data(x, y)
    particles.set_3d_properties(z)  
    return particles,

ani = animation.FuncAnimation(fig, animate, init_func=init, frames=num_steps, interval=20, blit=False)

plt.show()

ani.save('simul_3D.mp4', writer='ffmpeg', fps=100)

