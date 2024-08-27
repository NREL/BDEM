import numpy as np
import sys
    
lo = [-0.025, -0.025, 0.0]
hi = [0.025, 0.025, 0.1]

filename = "particle_input.dat"

radius = 0.001
mass = 430.0

buffer = 9.0 * radius

nx = 5
ny = 5
nz = 2

x_pos = np.linspace(lo[0] + 2.0 * buffer, hi[0] - 2.0 * buffer, nx)
y_pos = np.linspace(lo[1] + 2.0 * buffer, hi[1] - 2.0 * buffer, ny)
z_pos = np.linspace(lo[2] + 2.0 * buffer, hi[2] - 2.0 * buffer, nz)

z_pos, y_pos, x_pos = np.meshgrid(z_pos, y_pos, x_pos, indexing='ij')

pos_array = np.vstack((x_pos.flatten(), y_pos.flatten(), z_pos.flatten())).T

num_particles = np.shape(pos_array)[0]

with open(filename, "w") as fp:
    fp.write(f"     {num_particles:.0f}\n")
    
    for k, pos in enumerate(pos_array):
        fp.write("1 ")
        
        jitter = 0.0
        
        xx = pos[0] + jitter * (np.random.rand() - 0.5)
        yy = pos[1] + jitter * (np.random.rand() - 0.5)
        zz = pos[2] + jitter * (np.random.rand() - 0.5)
        
        fp.write(f"{xx:.6f} {yy:.6f} {zz:.6f} ") # x, y, z-position
        
        fp.write(f"{radius} {mass} ") # radius and mass
        
        fp.write(f"0.0 0.0 0.0 ") # Initial velocity
        
        elevation_angle = np.radians(90.0 - 5.0 * np.random.rand())
        azimuth_angle = np.radians(360.0 * np.random.rand())
        
        fp.write(f"0.0 {elevation_angle:.6f} {azimuth_angle:.6f} ") # angles of rotation
        
        fp.write(f"9\n") # Initial velocity

