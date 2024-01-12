import numpy as np
import matplotlib.pyplot as plt
import Non_LinearSolvers as solver

if __name__ == '__main__':

	domain = 4.0
	T = 2.0  
	nx = 41
	sigma = 1
	dx = domain/(nx-1)
	dt = sigma*dx
	nt = int(T/dt)
	c_coeff = 0.1
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
	u2 = solver.Richymyer(dt,dx,u_2,nt)
	u3 = solver.Lax_Friedrichs(dt,dx,u_3,nt)
	u4 = solver.MacCormack(dt,dx,u_4,nt)
	u5 = solver.Lax_Wendroff(dt,dx,u_5,nt)
	
	plt.plot(x,u1,label = 'Beam&Warming')
	plt.plot(x,u2,label = 'Richymyer')
	plt.plot(x,u3,label = 'Lax_Friedrichs')
	plt.plot(x,u4,label = 'MacCormack')
	plt.plot(x,u5,label = 'Lax_Wendroff')

	plt.xlabel('x(m)')
	plt.ylabel('u(m/s)')
	plt.ylim(None,1.7)
	plt.xlim(1,3.5)
	plt.title(f'Comparison with dx:{dx:.3f} , dt:{dt:.3f} with e = {c_coeff} and sigma: {sigma}')
	plt.legend()
	plt.show()
