import numpy as np
import matplotlib.pyplot as plt
import Non_LinearSolvers as solver

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
	
	u1 = solver.Beam_WarmingTA(nx,dt,dx,u_1,nt,True,c_coeff)
	"""
	u2 = solver.Beam_WarmingTA(nx,dt,dx,u_2,nt,True,0.025)
	u3 = solver.Beam_WarmingTA(nx,dt,dx,u_3,nt,True,0.05)
	u4 = solver.Beam_WarmingTA(nx,dt,dx,u_4,nt,True,0.075)
	u5 = solver.Beam_WarmingTA(nx,dt,dx,u_5,nt,True,0.1)
	u6 = solver.Beam_WarmingTA(nx,dt,dx,u_6,nt,True,0.125)
	"""
	
	#print(u)
	plt.plot(x,u1,label = 'u1')
	"""
	plt.plot(x,u2,label = 'u2')
	plt.plot(x,u3,label = 'u3')
	plt.plot(x,u4,label = 'u4')
	plt.plot(x,u5,label = 'u5')
	plt.plot(x,u6,label = 'u6')
	"""
	plt.xlabel('x(m)')
	plt.ylabel('u(m/s)')
	plt.ylim(-0.3,1.7)
	plt.xlim(1,3)
	plt.title(f'Comparison with dx:{dx} , dt:{dt} with e = {c_coeff}')
	plt.legend()
	plt.show()
