import Steps9_12 as sol

# Define Variables
nx = 50
ny = 50
nt  = 200
xmin = 0
xmax = 2
ymin = 0
ymax = 1

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)
# Initial Conditions

p  = sol.np.zeros((ny, nx))
b  = sol.np.zeros((ny, nx))
x  = sol.np.linspace(xmin, xmax, nx)
y  = sol.np.linspace(xmin, xmax, ny)

#Source
b[int(nx/4),int(ny/4)] = 100
b[int(3*nx/4),int(3*ny/4)] = -100

p = sol.Poisson2D(p,b,dx,dy,nt)

sol.plot2DLaplace(x,y,p)