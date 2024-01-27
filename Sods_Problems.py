import numpy as np
import NumericalSolvers_EulerEQN as sol
import matplotlib as plt
gamma = 1.4
nx = 20

P = np.zeros(nx)
u = np.zeros(nx)
rho = np.zeros(nx)

# Left side of domain

P[0:int(nx/2)] = 100000.0
u[0:int(nx/2)] = 0.0
rho[0:int(nx/2)] = 1.0

# Right side of Domain

P[int(nx/2):] = 10000.0
u[int(nx/2):] = 0.0
rho[int(nx/2):] = 0.125
"""
print('######################################################')
print(f'P matrix:\n{P}')
print('######################################################')
print(f'u matrix:\n{u}')
print('######################################################')
print(f'rho matrix:\n{rho}')
"""
et,q = sol.Conserved_VariableVec(rho,u,P,gamma)
"""
print('######################################################')
print(f'et matrix:\n{et}')
print('######################################################')
print(f'q matrix:\n{q}')
"""

F = sol.Flux_Vec(rho,u,P,et)

print(f'F matrix:\n{F}')

