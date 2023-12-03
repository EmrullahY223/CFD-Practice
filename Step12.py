import Steps9_12 as sol

##variable declarations
nx = 41
ny = 41
nt = 10
nit = 50 
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = sol.np.linspace(0, 2, nx)
y = sol.np.linspace(0, 2, ny)
X, Y = sol.np.meshgrid(x, y)


##physical variables
rho = 1
nu = .1
F = 1
dt = .01

#initial conditions
u = sol.np.zeros((ny, nx))

v = sol.np.zeros((ny, nx))

p = sol.np.ones((ny, nx))

b = sol.np.zeros((ny, nx))

u,v,p = sol.channel_flow_periodic(nt,u,v,dt,dx,dy,p,rho,nu,nx,ny,nit,F)

fig = sol.pyplot.figure()
sol.pyplot.contourf(X, Y, u, alpha=0.5, cmap=sol.cm.viridis)  
sol.pyplot.colorbar()
sol.pyplot.quiver(X[::3,::3],Y[::3,::3],u[::3,::3],v[::3,::3])
sol.pyplot.show()