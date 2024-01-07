import numpy as np
import matplotlib.pyplot as plt

def Upwind(nx,nt): # accurate for sigma <= 1 

    # define variables
    dx = 2/(nx-1)
    sigma = .1
    c = 1
    dt = sigma*dx/c

    u = np.zeros(nx)
    un = u.copy()

    x = np.linspace(0,2,nx)

    # Initial Conditions

    un[0] = 1
    u[0] = 1

    udiff = 1

    while udiff > 1e-5:
        un = u.copy()

        u[1:-1] = un[1:-1] - c*dt/dx*(un[1:-1]- un[0:-2])

        udiff = (np.sum(u) - np.sum(un))/np.sum(u)

    plt.plot(x,u)
    plt.show()

def LeapFrog(nx,nt): # accurate only at sigma = 1 !!!

    # define variables

    dx = 2/(nx-1)
    sigma = 1
    c = 1
    dt = sigma*dx/c

    u = np.zeros(nx)
    un_1 = u.copy()
    un = u.copy()

    # plotting aid 

    x = np.linspace(0,2,nx)

    # Initial Conditions

    un_1[0] = 1
    un[0] = 1
    u[0] = 1
    # Initial Calculation

    un[1:-1] = un_1[1:-1] - sigma*(un_1[1:-1]- un_1[0:-2])
    
    #print(un)
    print('*****************')

    # Leapfrog Calculations

    for n in range(nt+1):

        u[1:-1] = un_1[1:-1] - sigma*(un[2:]-un[0:-2])

        un_1 = un.copy()
        un = u.copy()

    print(u)
    plt.plot(x,u)
    plt.show()

def Lax_Friedrichs(nx,nt):#accurate only at sigma = 1 !!!
    # define variables
    dx = 2/(nx-1)
    sigma = 1
    c = 1
    dt = sigma*dx/c

    u = np.zeros(nx)
    un = u.copy()

    x = np.linspace(0,2,nx)

    # Initial Conditions

    un[0] = 1
    u[0] = 1

    for n in range(1,nt):
        un = u.copy()

        u[1:-1] = .5*(un[2:] + un[0:-2]) - c*dt/(2*dx)*(un[2:] - un[0:-2])

    plt.plot(x,u)
    plt.show()

def Lax_Wendroff(nx,nt):#accurate only at sigma = 1 !!!

    # define variables
    dx = 2/(nx-1)
    sigma = .001
    c = 2
    dt = sigma*dx/c
    print(f'dt : {dt}')

    u = np.zeros(nx)
    un = u.copy()

    x = np.linspace(0,2,nx)

    # Initial Conditions
    b = np.array(((c+c**2)/2 , 1-c**2 , (c**2-c)/2))
    un[0] = 1
    u[0] = 1

    for n in range(1,nt):
        un = u.copy()

        u[1:-1] = un[1:-1] - (((sigma)/2) * (un[2:] - un[0:-2])) + ((((sigma)**2)/2) * (un[2:] - 2*un[1:-1] + un[0:-2]))

        #u[1:-1] = b[0]*un[0:-2] + b[1]*un[1:-1] + b[2]*un[2:]

    print(u)
    plt.plot(x,u)
    plt.show()

if __name__ == '__main__':

    #LeapFrog(110,90)
    Upwind(31,100)
    #Lax_Friedrichs(100,90)
    #Lax_Wendroff(200,10000)
    print('hi')