"""
2D SPH written by Lalit Lal
"""

from math import sqrt, pi
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--mu', default=0.5)
parser.add_argument('--mass', default=1000)
parser.add_argument('--dist', default=20)
parser.add_argument('--scene', default='ss')
args = parser.parse_args()

# Simulation Parameters
N = 500  # Number of particles
width = 5 # x-axis
floor = 0. # y-axis
gravity = 9.81/1000 
max_vel = 2. 
wall_damp = 0.09 
floor_damp = 0.05
vel_damp = 0.5
dam = width # unless scene is specified as 'dam'

# particle stuff
mass = int(args.mass)
spacing = width/10  # starting space between particles
k = width/1000.  # gas constant
rest_d = 1.  # rest density
r = spacing * 1.75 # neighbor threshold
mu = float(args.mu)  # viscosity

# kernel stuff
h = int(args.dist)
poly6 = 315 / (64. * pi) / (h ** 9)
spiky_grad = -45 / pi / (h ** 6)
visc = 250
visc_lapl = visc * 45 / pi / (h ** 6)


def create_grid_particles(xmin, xmax, ymin, ymax, gap, count):
    result = []
    # box in centre of simulation window
    x, y = xmin, ymin
    for i in range(count):
        result.append([x,y])
        x += gap
        if x > xmax:
            x = xmin
            y += gap
        
    return result

def initialize():
    vel = []
    force = []
    pressure = []
    dens = []
    neighbor = []
    for i in range(N):
        vel.append([0., 0.])
        force.append([0., 0.])
        pressure.append(0.)
        dens.append(0.)
        neighbor.append([]) # start off with empty list of neighbors for each particle

    # position is already initialized from particle creation
    return vel, force, pressure, dens, neighbor

initial = []
scene = args.scene
# sinlgle stream
if scene == 'ss':
    initial = create_grid_particles(-width/4, width/4, width/4, width/4, width/10, N)
# multistream
elif scene == 'ms':
    initial = create_grid_particles(-width/2, -width/4, width/4, width/4, width/10, int(N/2))
    initial += create_grid_particles(width/4, width/2, width/4, width/4, width/10, int(N/2))
# dam break
elif scene == 'dam':
    initial = create_grid_particles(-width, -width/4, width/4, width/4, width/10, N)
    dam = width/8.
#default
else:
    initial = create_grid_particles(-width/4, width/4, width/4, width/4, width/10, N)


# Initial dam break conditions
curr_pos = np.asarray(initial)
prev_pos = np.asarray(initial)
visual = np.asarray(initial)
velocity, force, pressure, dens, neighbor = initialize()

def pressure_projection():
    # first build pressure from densities
    for i in range(N):
        pressure[i] = k * (dens[i] - rest_d)

    # force due to pressure
    for i in range(N):
        dPress = [0, 0]
        for n in neighbor[i]:
            ni, q = n[0], n[1]
            rij = [curr_pos[ni][0] - curr_pos[i][0], curr_pos[ni][1] - curr_pos[i][1]]
            rij_mag = sqrt((rij[0]) ** 2 + (rij[1]) ** 2)
            rij_normd = [rij[0]/rij_mag, rij[1]/rij_mag]
        
            p = [0,0]
            p[0] = rij_normd[0] * mass * (pressure[i] + pressure[ni]) / (2 * dens[ni]) * spiky_grad * (h - rij_mag)**2
            p[1] = rij_normd[1] * mass * (pressure[i] + pressure[ni]) / (2 * dens[ni]) * spiky_grad * (h - rij_mag)**2
            P = [rij_normd[0] * p[0], rij_normd[1] * p[1]]

            dPress = [dPress[0] + P[0], dPress[1] + P[1]]
            force[ni][0] += P[0]
            force[ni][1] += P[1]
        force[i][0] -= dPress[0]
        force[i][1] -= dPress[1]

def build_neighbors_and_density():
    for i in range(N):
        d = 0.
        for j in range(N):
            if (i < j):
                dist = sqrt((curr_pos[i][0] - curr_pos[j][0]) ** 2 + (curr_pos[i][1] - curr_pos[j][1]) ** 2)
                sq_dist = dist**2
                if (sq_dist < r**2):
                    dj = mass * poly6 * (h**2 - sq_dist) ** 3
                    d += dj
                    dens[j] += dj
                    neighbor[i].append([j,dj])
        
        dens[i] += d

def apply_viscosity_force():
    for i in range(N):
        for j, n in enumerate(neighbor[i]):
            ni, _ = n[0], n[1]
            rij = [curr_pos[ni][0] - curr_pos[i][0], curr_pos[ni][1] - curr_pos[i][1]]
            l = sqrt(rij[0]**2 + rij[1]**2)
            q = l / r
            rij_n = [rij[0]/l, rij[1]/l]
            v_diff = [velocity[i][0] - velocity[ni][0], velocity[i][1] - velocity[ni][1]]
            u = v_diff[0] * rij_n[0] + v_diff[1] * rij_n[1]
            
            if (u > 0):
                base = (1 - q) * (mu * u)
                Fv = [base * rij_n[0], base * rij_n[1]]
                velocity[i][0] -= Fv[0] * 0.5
                velocity[i][1] -= Fv[1] * 0.5
                velocity[ni][0] += Fv[0] * 0.5
                velocity[ni][1] += Fv[1] * 0.5

def advect_particles():
    for i in range(N):
        curr_pos[i][0] += velocity[i][0]
        curr_pos[i][1] += velocity[i][1]

def apply_body_forces():
    for i in range(N):
        curr_pos[i][0] += force[i][0]
        curr_pos[i][1] += force[i][1]
        visual[i] = curr_pos[i]
        force[i] = [0, -gravity]
    
def save_state():
    for i in range(N):
        prev_pos[i] = curr_pos[i]

def apply_velocity_boundary_conditions():
    for i in range(N):
        velocity[i][0] = curr_pos[i][0] - prev_pos[i][0]
        velocity[i][1] = curr_pos[i][1] - prev_pos[i][1]
        V_mag = sqrt(velocity[i][0] ** 2 + velocity[i][1] ** 2)
        
        if V_mag > max_vel:
            velocity[i][0] *= vel_damp
            velocity[i][1] *= vel_damp
        
        if curr_pos[i][0] < -width:
            force[i][0] -= (curr_pos[i][0] - -width) * wall_damp
            visual[i][0] = -width
        
        if curr_pos[i][0] > dam:
            force[i][0] -= (curr_pos[i][0] - dam) * wall_damp
        
        if curr_pos[i][0] > dam:
            visual[i][0] = dam
        
        if curr_pos[i][1] < floor:
            force[i][1] -= (curr_pos[i][1] - width) * floor_damp # floor might need to change
            visual[i][1] = floor


def step():
    global curr_pos, prev_pos, velocity, force, pressure, dens, neighbor

    save_state()
    advect_particles()
    apply_body_forces()
    apply_velocity_boundary_conditions()    
    # clear state
    for i in range(N):
        dens[i] = 0 
        neighbor[i] = []

    # # # DENSITY # # #
    build_neighbors_and_density()
    # # # END DENSITY

    # # # PRESSURE
    pressure_projection()
    # # # END PRESSURE

    # # # VISCOSITY
    apply_viscosity_force()
    # # # END VISCOSITY


# Build Simulation View
fig = plt.figure()
axes = fig.add_subplot(xlim=(-width, width), ylim=(floor, 2 * width))
points, = axes.plot([], [], 'bo', ms=20)
frame = 0

def animate(i):
    global curr_pos, dam, frame
    step()
    frame += 1
    if frame == 150:  # Break the dam at frame 150
        dam = width
    points.set_data(visual[:, 0], visual[:, 1])  # Position updates of actual visual
    return points,

mywriter = animation.FFMpegWriter(fps=30)
ani = animation.FuncAnimation(fig, animate, interval=10, blit=True, save_count = 250)
filename = '{visc}_{scene}.'.format(visc=mu, scene=scene) + 'mp4'
ani.save(filename, writer=mywriter)
print('Video Saved!')

