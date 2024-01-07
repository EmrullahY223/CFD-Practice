import Steps9_12 as sol

# Define Variables
nx = 51
ny = 51
c = 1
dx = 2/(nx-1)
dy = 1/(ny-1)

# Initial Conditions

p = sol.np.zeros((ny,nx))

#Plotting Aids

x = sol.np.linspace(0,2,nx)
y = sol.np.linspace(0,1,ny)

# Boundary Conditions

p[:, 0] = 0  # p = 0 @ x = 0
p[:, -1] = y  # p = y @ x = 2
p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

#sol.plot2DLaplace(x,y,p)

p = sol.laplace2d(p,y,dx,dy,1e-6)

sol.plot2DLaplace(x,y,p)