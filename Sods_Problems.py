import numpy as np
import NumericalSolvers_EulerEQN as sol
import matplotlib.pyplot as plt

# Initial Conditions #

# Right Side

rhoR = 0.125
uR = 0.0
PR = 10.0e3

# Left Side

rhoL = 1.0
uL = 0.0
PL = 100.0e3

# Solution Settings #

xmin = -10.0; xmax = 10.0; Tfin = 0.01

nx = 50

dx = (xmax - xmin)/nx

CFL = 0.4; gamma = 1.4
a = np.sqrt(gamma * PL/rhoL)

dt = CFL*(dx/a)

print('##############################################')
print(f'dt = {dt}')
print('##############################################')
print(f'dx = {dx}')
print('##############################################')
print(f'dt/dx = {dt/dx}')
print('##############################################')
print(f'a = {a}')

# Initial Arrays

P = np.zeros(nx)
u = np.zeros(nx)
rho = np.zeros(nx)

# Array for plotting

x = np.linspace(-10,10,nx)

# Left side of Domain

P[0:int(nx/2)] = PL
u[0:int(nx/2)] = uL
rho[0:int(nx/2)] = rhoL

# Right side of Domain

P[int(nx/2):] = PR
u[int(nx/2):] = uR
rho[int(nx/2):] = rhoR
#####################################################
# Testing

# q = sol.Conserved_VariableVec(rho,u,P,gamma)

# F = sol.Flux_Vec(q,gamma)

####################################################

# Solution and Post-Processing

q0 = sol.Conserved_VariableVec(rho,u,P,gamma)

Density, Velocity, Pressure, SpeedOfSound, MachNumber = sol.RichtMyer(q0,dt,dx,Tfin,gamma)

print('##############################################')
print(f'Density = \n{Density}')
print('##############################################')
print(f'velocity = \n{Velocity}')
print('##############################################')
print(f'Pressure = \n{Pressure}')
print('##############################################')
print(f'SpeedOfSound = \n{SpeedOfSound}')
print('##############################################')
print(f'MachNumber = \n{MachNumber}')

# Plot The Results

fig, axs = plt.subplots(2,3)
axs[0,0].plot(x,Pressure)
axs[0,0].set_title('Pressure')
axs[0,0].set(xlabel = 'x(m)', ylabel = 'kN/(m^2)')

axs[0,1].plot(x,Velocity)
axs[0,1].set_title('Velocity')
axs[0,1].set(xlabel = 'x(m)', ylabel = 'm/s')

axs[0,2].plot(x,SpeedOfSound)
axs[0,2].set_title('Speed of Sound')
axs[0,2].set(xlabel = 'x(m)', ylabel = 'm/s')

axs[1,0].plot(x,Density)
axs[1,0].set_title('Density')
axs[1,0].set(xlabel = 'x(m)', ylabel = 'kg/(m^3)')

axs[1,1].plot(x,MachNumber)
axs[1,1].set_title('Mach Number')
axs[1,1].set(xlabel = 'x(m)', ylabel = 'Ma')

plt.show()