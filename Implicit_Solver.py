import numpy as np
import matplotlib.pyplot as plt

def Beam_WarmingTA(nx,dt,dx,u,nt):

	# Initialization of arrays

	a = np.zeros(nx)
	b = np.ones(nx)		#b matrix is 1 for this problem
	c = np.zeros(nx)
	c_star = np.zeros(nx)
	d = np.zeros(nx)		#d matrix is un for this problem
	d_star = np.zeros(nx)
	un = np.zeros(nx)

	# Calculation of Velocity Field
	
	for n in range(nt):
		
		un[:] = u[:]

		a[0] = -(dt/(4*dx)) * un[0]
		a[1:] = -(dt/(4*dx)) * un[0:-1]

		c[:-1] = (dt/(4*dx)) * un[1:]

		c_star[0] = c[0]/b[0]
		c_star[1:-1] = c[1:-1]/(b[1:-1] - a[1:-1]*c_star[0:-2])

		d[:] = un[:]
		d[0] = d[0] - a[0]*un[0]

		d_star[0] = d[0]/b[0]
		d_star[1:] = (d[1:] - a[1:]*d_star[0:-1])/(b[1:] - a[1:]*c_star[0:-1])

		u[nx-1] = d_star[nx-1]

		for i in range(nx-2, -1, -1):
			
			u[i] = d_star[i] - c_star[i] * un[i+1]

		#print(u)
		u[0] = 1
	"""
	
	for n in range(nt):
		un[:] = u[:]

		a[0] = -(dt/(4*dx)) * un[0]

		for i in range(1,nx):
			a[i] = -dt/(4*dx)*un[i-1]
		for i in range(nx-1):
			c[i] = dt/(dx*4)*un[i+1]
		for i in range(nx):
			d[i] = un[i]

		d[0] = d[0] - a[0]*un[0]
		c_star[0] = c[0]/b[0]
		d_star[0] = d[0]/b[0]

		for i in range(1,nx-1):
			c_star[i] = c[i]/(b[i] - c_star[i-1]*a[i])

		for i in range(1,nx):
			d_star[i] = (d[i] - d_star[i-1]*a[i])/(b[i] - c_star[i-1]*a[i])
		
		u[nx-1] = d_star[nx-1]

		for i in range(nx-2, -1, -1):
			
			u[i] = d_star[i] - c_star[i] * un[i+1]

		print(u)
		u[0] = 1
	"""
	return u

def Damping(c_coeff, un):
	
	D = -c_coeff*(un[4:] - 4*un[3:-1] + 6*un[2:-2] - 4*un[1:-3] + un[0:-4])

	return D

if __name__ == '__main__':

	domain = 4.0
	T = 0.4  
	nx = 41
	sigma = 1.0
	c_coeff = 0.0
	dx = domain/(nx-1)
	dt = sigma*dx
	nt = int(T/dt)
	print(f'dx: {dx} , dt: {dt}')

	#IC's and BC's

	u = np.zeros(nx)
	x = np.linspace(0,domain,nx)

	for i in range(nx):

		if x[i] <= 2:

			u[i] = 1

		else:

			u[i] = 0

	# Solution Iterations
	#print(u)

	u = Beam_WarmingTA(nx,dt,dx,u,nt)

	#print(u)
	plt.plot(x,u)
	plt.show()
