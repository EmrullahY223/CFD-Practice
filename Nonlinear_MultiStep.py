import numpy as np

def Initial_and_boudary_conditions(nx,nt):
    u = np.zeros((nt,nx))

    # Boundary Conditions

    u[:,0] = 1.0
    u[:,-1] = 0.0

    # Initial Conditions

    half = int(nx/2)

    u[0,:half] = 1.0

    return u

def analytical(nx, tmax, xmax):

    dx = xmax/(nx-1)
    dt = dx/0.5
    nt = int(tmax/dt+1)

    # Initial and boundary conditions

    u = Initial_and_boudary_conditions(nx,nt)

    x = np.linspace(0,xmax,nx)

    for n in range(nt-1):

        u[n+1, 1:-1] = u[n, 0:-2]

    return u,x,nt

def plot (u,x,nt,u_analitical,x_analytical,nt_analytical,step,dt,dt_analytical):

    import matplotlib.pyplot as plt
    
    plt.figure()
    ax = plt.subplot(1,1,1)

    for n in range(0,nt,step):
        t = n*dt
        ax.plot(x,u[n,:],linestyle = '-', label = 't:' + str(t))

    t_analytical = dt_analytical* (nt_analytical-1)

    ax.plot(x_analytical,u_analitical[-1,:],linestyle = '--',label = 't:' + str(t_analytical))
    ax.legend( bbox_to_anchor=(1.02,1), loc=2)
    plt.xlabel('x (m)')
    plt.ylabel('u (m/s)')
    plt.show()

def plotCFL (u1,u2,u3,x,u_anal,x_anal,nt1,nt2,nt3,nt_anal,dt1,dt2,dt3,dt_anal,sigma1,sigma2,sigma3):

    import matplotlib.pyplot as plt
    
    plt.figure()
    ax = plt.subplot(1,1,1)

    ax.plot(x,u1[nt1-1,:],linestyle = '-', label = 't:' + str(2) + ' sigma:'+str(sigma1))
    ax.plot(x,u2[nt2-1,:],linestyle = '-', label = 't:' + str(2) + ' sigma:'+str(sigma2))
    ax.plot(x,u3[nt3-1,:],linestyle = '-', label = 't:' + str(2) + ' sigma:'+str(sigma3))

    t_analytical = dt_anal* (nt_anal-1)

    ax.plot(x_anal,u_anal[-1,:],linestyle = '--',label = 't:' + str(t_analytical))
    ax.legend( bbox_to_anchor=(1.02,1), loc=2)
    plt.xlabel('x (m)')
    plt.ylabel('u (m/s)')
    plt.show()

def plotComp(u1,u2,u3,x,u_anal,x_anal,nt,nt_anal,dt,dt_anal,sigma):
    
    import matplotlib.pyplot as plt
    
    plt.figure()
    ax = plt.subplot(1,1,1)
    
    ax.plot(x,u1[nt-1,:], linestyle = '-', label = 'Scheme: RichtMyer' + ' sigma:'+str(sigma))
    ax.plot(x,u2[nt-1,:], linestyle = '-', label = 'Scheme: Lax-Friedrichs' + ' sigma:'+str(sigma))
    ax.plot(x,u3[nt-1,:], linestyle = '-', label = 'Scheme: MacCormack' + ' sigma:'+str(sigma))
    
    ax.plot(x_anal,u_anal[-1,:],linestyle = '--',label = 'AnaLytical')
    ax.legend( bbox_to_anchor=(1.02,1), loc=2)
    plt.xlabel('x (m)')
    plt.ylabel('u (m/s)')
    plt.show()

def RichtMyer(sigma,nx,xmax,tmax,step):

    dx = xmax/(nx-1)
    dt = sigma*dx
    nt = int(tmax/dt+1)
    print()

    # IC's and BC's

    u = Initial_and_boudary_conditions(nx,nt)

    # Plot aids
    x = np.linspace(0,xmax,nx)

    # Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    for n in range(nt-1):

        uh = u.copy()

        uh[n, 1:-1] = 0.5 * (u[n,2:] + u[n,0:-2]) - dt/(4*dx) * (F(u[n,2:]) - F(u[n,0:-2]))
        #print(uh)

        u[n+1, 1:-1] = u[n,1:-1] - dt/(dx*2) * (F(uh[n,2:]) - F(uh[n,0:-2]))

    u_anal, x_anal, nt_anal = analytical(nx,tmax,xmax)

    dt_anal = 2*dx
    
    return u,x,nt,dt,dt_anal

    #plot(u,x,nt,u_anal,x_anal,nt_anal,step,dt,dt_anal)

def RichtMyer_cfl(sigma1,sigma2,sigma3,nx,xmax,tmax):

    dx = xmax/(nx-1)
    dt1 = sigma1*dx
    dt2 = sigma2*dx
    dt3 = sigma3*dx
    nt1 = int((tmax/dt1)+1)
    nt2 = int((tmax/dt2)+1)
    nt3 = int((tmax/dt3)+1)
    print(f'nt1: {nt1}')
    print(f'nt2: {nt2}')
    print(f'nt3: {nt3}')

    # IC's and BC's

    u1 = Initial_and_boudary_conditions(nx,nt1)
    u2 = Initial_and_boudary_conditions(nx,nt2)
    u3 = Initial_and_boudary_conditions(nx,nt3)

    # Plot aids
    x = np.linspace(0,xmax,nx)

    # Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations

    # Solution Iteration for Sigma 1
    for n in range(nt1-1):

        uh = u1.copy()

        uh[n, 1:-1] = 0.5 * (u1[n,2:] + u1[n,0:-2]) - dt1/(4*dx) * (F(u1[n,2:]) - F(u1[n,0:-2]))

        u1[n+1, 1:-1] = u1[n,1:-1] - dt1/(dx*2) * (F(uh[n,2:]) - F(uh[n,0:-2]))
    
    # Solution Iteration for Sigma 2 
    for n in range(nt2-1):

        uh = u2.copy()

        uh[n, 1:-1] = 0.5 * (u2[n,2:] + u2[n,0:-2]) - dt2/(4*dx) * (F(u2[n,2:]) - F(u2[n,0:-2]))

        u2[n+1, 1:-1] = u2[n,1:-1] - dt2/(dx*2) * (F(uh[n,2:]) - F(uh[n,0:-2]))

    # Solution Iteration for Sigma 3
    for n in range(nt3-1):

        uh = u3.copy()

        uh[n, 1:-1] = 0.5 * (u3[n,2:] + u3[n,0:-2]) - dt3/(4*dx) * (F(u3[n,2:]) - F(u3[n,0:-2]))

        u3[n+1, 1:-1] = u3[n,1:-1] - dt3/(dx*2) * (F(uh[n,2:]) - F(uh[n,0:-2]))

    u_anal, x_anal, nt_anal = analytical(nx,tmax,xmax)

    dt_anal = 2*dx

    plotCFL(u1,u2,u3,x,u_anal,x_anal,nt1,nt2,nt3,nt_anal,dt1,dt2,dt3,dt_anal,sigma1,sigma2,sigma3)

def Lax_Friedrichs_MultiStep(sigma1,nx,xmax,tmax,step):
    
    # Declaring Variables
    
    dx = xmax/(nx-1)
    dt1 = sigma1*dx
    nt1 = int((tmax/dt1)+1)
    print(f'nt1: {nt1}')

    # IC's and BC's

    u1 = Initial_and_boudary_conditions(nx,nt1)

    # Plot aids
    x = np.linspace(0,xmax,nx)

    # Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    
    # Solution Iteration for Sigma 1
    
    for n in range(nt1 - 1):
        
        unh_p = u1.copy()
        unh_n = u1.copy()
        
        unh_p[n,1:-1] = 0.5*(u1[n,1:-1] + u1[n,2:]) - dt1/(dx*2)*(F(u1[n,2:]) - F(u1[n,1:-1]))
        unh_n[n,1:-1] = 0.5*(u1[n,1:-1] + u1[n,0:-2]) - dt1/(dx*2)*(F(u1[n,1:-1]) - F(u1[n,0:-2]))
        
        u1[n+1,1:-1] = u1[n,1:-1] - dt1/dx*(F(unh_p[n,1:-1]) - F(unh_n[n,1:-1]))
        
    
    u_anal, x_anal, nt_anal = analytical(nx,tmax,xmax)

    dt_anal = 2*dx
    
    return u1

    #plot(u1,x,nt1,u_anal,x_anal,nt_anal,step,dt1,dt_anal)

def Lax_Friedrichs_MultiStepCFL(sigma1,sigma2,sigma3,nx,xmax,tmax):
    
    # Declaring Variables
    
    dx = xmax/(nx-1)
    dt1 = sigma1*dx
    dt2 = sigma2*dx
    dt3 = sigma3*dx
    nt1 = int((tmax/dt1)+1)
    nt2 = int((tmax/dt2)+1)
    nt3 = int((tmax/dt3)+1)
    print(f'nt1: {nt1}')
    print(f'nt2: {nt2}')
    print(f'nt3: {nt3}')

    # IC's and BC's

    u1 = Initial_and_boudary_conditions(nx,nt1)
    u2 = Initial_and_boudary_conditions(nx,nt2)
    u3 = Initial_and_boudary_conditions(nx,nt3)

    # Plot aids
    x = np.linspace(0,xmax,nx)

    # Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    
    # Solution Iteration for Sigma 1
    
    for n in range(nt1 - 1):
        
        unh_p = u1.copy()
        unh_n = u1.copy()
        
        unh_p[n,1:-1] = 0.5*(u1[n,1:-1] + u1[n,2:]) - dt1/(dx*2)*(F(u1[n,2:]) - F(u1[n,1:-1]))
        unh_n[n,1:-1] = 0.5*(u1[n,1:-1] + u1[n,0:-2]) - dt1/(dx*2)*(F(u1[n,1:-1]) - F(u1[n,0:-2]))
        
        u1[n+1,1:-1] = u1[n,1:-1] - dt1/dx*(F(unh_p[n,1:-1]) - F(unh_n[n,1:-1]))
        
    # Solution Iteration for Sigma 2
    
    for n in range(nt2 - 1):
        
        unh_p = u2.copy()
        unh_n = u2.copy()
        
        unh_p[n,1:-1] = 0.5*(u2[n,1:-1] + u2[n,2:]) - dt2/(dx*2)*(F(u2[n,2:]) - F(u2[n,1:-1]))
        unh_n[n,1:-1] = 0.5*(u2[n,1:-1] + u2[n,0:-2]) - dt2/(dx*2)*(F(u2[n,1:-1]) - F(u2[n,0:-2]))
        
        u2[n+1,1:-1] = u2[n,1:-1] - dt2/dx*(F(unh_p[n,1:-1]) - F(unh_n[n,1:-1]))
        
    # Solution Iteration for Sigma 3
    
    for n in range(nt3 - 1):
        
        unh_p = u3.copy()
        unh_n = u3.copy()
        
        unh_p[n,1:-1] = 0.5*(u3[n,1:-1] + u3[n,2:]) - dt3/(dx*2)*(F(u3[n,2:]) - F(u3[n,1:-1]))
        unh_n[n,1:-1] = 0.5*(u3[n,1:-1] + u3[n,0:-2]) - dt3/(dx*2)*(F(u3[n,1:-1]) - F(u3[n,0:-2]))
        
        u3[n+1,1:-1] = u3[n,1:-1] - dt3/dx*(F(unh_p[n,1:-1]) - F(unh_n[n,1:-1]))
        
    u_anal, x_anal, nt_anal = analytical(nx,tmax,xmax)

    dt_anal = 2*dx

    plotCFL(u1,u2,u3,x,u_anal,x_anal,nt1,nt2,nt3,nt_anal,dt1,dt2,dt3,dt_anal,sigma1,sigma2,sigma3)

def MacCormack(sigma,nx,xmax,tmax,step):
    
    dx = xmax/(nx-1)
    dt = sigma*dx
    nt = int(tmax/dt+1)
    print()

    # IC's and BC's

    u = Initial_and_boudary_conditions(nx,nt)

    # Plot aids
    x = np.linspace(0,xmax,nx)

    # Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    
    for n in range(nt-1):
        
        us = u.copy()
        
        us[n,1:-1] = u[n,1:-1] - dt/dx*(F(u[n,2:]) - F(u[n,1:-1]))
        
        u[n+1,1:-1] = 0.5* ((u[n,1:-1] + us[n,1:-1]) - dt/dx*(F(us[n,1:-1]) - F(us[n,0:-2])))
        
    u_anal, x_anal, nt_anal = analytical(nx,tmax,xmax)

    dt_anal = 2*dx
    
    return u
    
    #plot(u,x,nt,u_anal,x_anal,nt_anal,step,dt,dt_anal)
    
def MacCormackCFL(sigma1,sigma2,sigma3,nx,xmax,tmax):
    
    dx = xmax/(nx-1)
    dt1 = sigma1*dx
    dt2 = sigma2*dx
    dt3 = sigma3*dx
    nt1 = int((tmax/dt1)+1)
    nt2 = int((tmax/dt2)+1)
    nt3 = int((tmax/dt3)+1)
    print(f'nt1: {nt1}')
    print(f'nt2: {nt2}')
    print(f'nt3: {nt3}')

    # IC's and BC's

    u1 = Initial_and_boudary_conditions(nx,nt1)
    u2 = Initial_and_boudary_conditions(nx,nt2)
    u3 = Initial_and_boudary_conditions(nx,nt3)

    # Plot aids
    x = np.linspace(0,xmax,nx)

    # Flux fucn

    F = lambda u: 0.5*(u**2)
    
    # Solution Iterations
    
    # Solution Iteration for Sigma 1
    
    for n in range(nt1-1):
        
        us = u1.copy()
        
        us[n,1:-1] = u1[n,1:-1] - dt1/dx*(F(u1[n,2:]) - F(u1[n,1:-1]))
        
        u1[n+1,1:-1] = 0.5* ((u1[n,1:-1] + us[n,1:-1]) - dt1/dx*(F(us[n,1:-1]) - F(us[n,0:-2])))
        
    # Solution Iteration for Sigma 2
    
    for n in range(nt2-1):
        
        us = u2.copy()
        
        us[n,1:-1] = u2[n,1:-1] - dt2/dx*(F(u2[n,2:]) - F(u2[n,1:-1]))
        
        u2[n+1,1:-1] = 0.5* ((u2[n,1:-1] + us[n,1:-1]) - dt2/dx*(F(us[n,1:-1]) - F(us[n,0:-2])))
        
    # Solution Iteration for Sigma 3
    
    for n in range(nt3-1):
        
        us = u3.copy()
        
        us[n,1:-1] = u3[n,1:-1] - dt3/dx*(F(u3[n,2:]) - F(u3[n,1:-1]))
        
        u3[n+1,1:-1] = 0.5* ((u3[n,1:-1] + us[n,1:-1]) - dt3/dx*(F(us[n,1:-1]) - F(us[n,0:-2])))
       
    # Analytical Solution
    
    u_anal, x_anal, nt_anal = analytical(nx,tmax,xmax)

    dt_anal = 2*dx

    plotCFL(u1,u2,u3,x,u_anal,x_anal,nt1,nt2,nt3,nt_anal,dt1,dt2,dt3,dt_anal,sigma1,sigma2,sigma3)
    
def Comparison(sigma,nx,xmax,tmax,step):
    
    u1,x,nt,dt,dt_anal = RichtMyer(sigma,nx,xmax,tmax,step)
    u2 = Lax_Friedrichs_MultiStep(sigma,nx,xmax,tmax,step)
    u3 = MacCormack(sigma,nx,xmax,tmax,step)
    u_anal, x_anal, nt_anal = analytical(nx,tmax,xmax)
    
    plotComp(u1, u2, u3, x, u_anal, x_anal, nt, nt_anal, dt, dt_anal, sigma)

if __name__ == '__main__':

    #RichtMyer(0.25,101,4.0,2.0,50)
    #RichtMyer_cfl(sigma1 = 2.0,sigma2 = 1.5,sigma3 = 1.0,nx = 101,xmax = 10.0,tmax = 5.0)
    #Lax_Friedrichs_MultiStepCFL(sigma1 = 1.,sigma2 = 1,sigma3 = .5,nx = 101,xmax = 4.0,tmax = 2.0)
    #Lax_Friedrichs_MultiStep(sigma1 = 0.5,nx = 101,xmax = 4.0,tmax = 2.0,step = 24)
    #MacCormack(0.5, 101, 4.0, 2.0, 24)
    #MacCormackCFL(sigma1 = 1,sigma2 = 0.5,sigma3 = 0.25,nx = 101,xmax = 4.0,tmax = 2.0)
    Comparison(.5, 101, 4.0, 2.0, 1)