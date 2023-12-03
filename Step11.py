import Steps9_12 as sol

nx = 41
ny = 41
nt = 700
nit = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = sol.np.linspace(0, 2, nx)
y = sol.np.linspace(0, 2, ny)
X, Y = sol.np.meshgrid(x, y)

rho = 1
nu = .1
dt = .001

u = sol.np.zeros((ny, nx))
v = sol.np.zeros((ny, nx))
p = sol.np.zeros((ny, nx)) 
b = sol.np.zeros((ny, nx))

u,v,p = sol.cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu,nx,ny,nit)

fig = sol.pyplot.figure(figsize=(11,7), dpi=100)
# plotting the pressure field as a contour
sol.pyplot.contourf(X, Y, p, alpha=0.5, cmap=sol.cm.viridis)  
sol.pyplot.colorbar()
# plotting the pressure field outlines
sol.pyplot.contour(X, Y, p, cmap=sol.cm.viridis)  
# plotting velocity field
sol.pyplot.streamplot(X, Y, u, v) 
sol.pyplot.xlabel('X')
sol.pyplot.ylabel('Y')

sol.pyplot.show()
