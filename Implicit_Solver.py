import numpy as np
import matplotlib.pyplot as plt

def Beam_WarmingTA(nx,dt,dx,u,nt,damping,c_coeff):

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
		for i in range(1,nx-1):
			c_star[i] = c[i]/(b[i] - (a[i]*c_star[i-1]))
		
		# Introduce damping if needed
			
		if(damping):
			d[0] = un[0]
			d[1] = un[1]

			d[2:-2] = un[2:-2] + Damping(c_coeff, un)

			d[nx-2] = un[nx-2]
			d[nx-1] = un[nx-1]

		else:
			d[:] = un[:]
		
		d[0] = d[0] - a[0]*un[0]
		# *********************************
		print(f'd: {d}')	
		d_star[0] = d[0]/b[0]
		for i in range(1,nx):
			d_star[i] = (d[i] - a[i]*d_star[i-1])/(b[i] - a[i]*c_star[i-1])
		
		u[nx-1] = d_star[nx-1]

		for i in range(nx-2, -1, -1):
			
			u[i] = d_star[i] - c_star[i] * un[i+1]

		#print(f"u = {u}")
		u[0] = 1
	
	return u

def Richymyer(dt,dx,u,nt):
	# Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    for n in range(nt):

        uh = u.copy()

        uh[1:-1] = 0.5 * (u[2:] + u[0:-2]) - dt/(4*dx) * (F(u[2:]) - F(u[0:-2]))
        #print(uh)

        u[1:-1] = u[1:-1] - dt/(dx*2) * (F(uh[2:]) - F(uh[0:-2]))
	

def Lax_Friedrichs(dt,dx,u,nt):

	# Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    
    for n in range(nt):
        
        unh_p = u.copy()
        unh_n = u.copy()
        
        unh_p[1:-1] = 0.5*(u[1:-1] + u[2:]) - dt/(dx*2)*(F(u[2:]) - F(u[1:-1]))
        unh_n[1:-1] = 0.5*(u[1:-1] + u[0:-2]) - dt/(dx*2)*(F(u[1:-1]) - F(u[0:-2]))
        
        u[1:-1] = u[1:-1] - dt/dx*(F(unh_p[1:-1]) - F(unh_n[1:-1]))
    
	return u

def MacCormack(dt,dx,u,nt):

	# Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    
    for n in range(nt):
        
        us = u.copy()
        
        us[1:-1] = u[1:-1] - dt/dx*(F(u[2:]) - F(u[1:-1]))
        
        u[1:-1] = 0.5* ((u[1:-1] + us[1:-1]) - dt/dx*(F(us[1:-1]) - F(us[0:-2])))

	return u

def Lax_Wendroff(dt,dx,u,nt):

	# Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations

	for n in range(nt):

		un = u.copy()

		u[1:-1] = un[1:-1] - dt/(2*dx)*(F(un[2:] - F(u[0:-2]))) + ...
		(dt**2)/(4*dx**2)*((un[2:]+un[1:-1])*(F(un[2:])-F(un[1:-2])) - (un[1:-1] + un[0:-2])*(F(un[1:-1])-F(un[0:-2])))

	return u

def Damping(c_coeff, un):
	
	D = -c_coeff*(un[4:] - 4*un[3:-1] + 6*un[2:-2] - 4*un[1:-3] + un[0:-4])

	return D

if __name__ == '__main__':

	domain = 4.0
	T = 0.4  
	nx = 41
	sigma = 1.0
	dx = domain/(nx-1)
	dt = sigma*dx
	nt = int(T/dt)
	c_coeff = 0.125
	print(f'dx: {dx} , dt: {dt}')

	#IC's and BC's

	u_1 = np.zeros(nx)
	u_2 = np.zeros(nx)
	u_3 = np.zeros(nx)
	u_4 = np.zeros(nx)
	u_5 = np.zeros(nx)
	u_6 = np.zeros(nx)
	x = np.linspace(0,domain,nx)

	for i in range(nx):

		if x[i] <= 2:

			u_1[i] = 1
			u_2[i] = 1
			u_3[i] = 1
			u_4[i] = 1
			u_5[i] = 1
			u_6[i] = 1
		else:

			u_1[i] = 0
			u_2[i] = 0
			u_3[i] = 0
			u_4[i] = 0
			u_5[i] = 0
			u_6[i] = 0

	# Solution Iterations
	
	u1 = Beam_WarmingTA(nx,dt,dx,u_1,nt,False,0)
	u2 = Beam_WarmingTA(nx,dt,dx,u_1,nt,True,0.25)
	u3 = Beam_WarmingTA(nx,dt,dx,u_1,nt,True,0.5)
	u4 = Beam_WarmingTA(nx,dt,dx,u_1,nt,True,0.75)
	u5 = Beam_WarmingTA(nx,dt,dx,u_1,nt,True,0.1)
	u6 = Beam_WarmingTA(nx,dt,dx,u_1,nt,True,0.125)
	#print(u)
	plt.plot(x,u1,label = 'u1')
	plt.plot(x,u2,label = 'u2')
	plt.plot(x,u3,label = 'u3')
	plt.plot(x,u4,label = 'u4')
	plt.plot(x,u5,label = 'u5')
	plt.plot(x,u6,label = 'u6')
	plt.xlabel('x(m)')
	plt.ylabel('u(m/s)')
	plt.legend()
	plt.show()
