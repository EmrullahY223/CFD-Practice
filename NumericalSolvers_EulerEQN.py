import numpy as np

def Conserved_VariableVec(rho,u,P,gamma):

    et = P/((gamma-1)*rho)+ 1/2*u**2
    
    q = np.array([rho,rho*u,rho*et])

    return q

def Flux_Vec(q,gamma):
    
    rho = q[0]
    u = q[1]/rho
    et = q[2]/rho
    P = (gamma-1)*rho*(et-0.5*u**2)

    F1 = rho*u
    F2 = (rho*u**2) + P
    F3 = (rho*et + P) * u

    F = np.array([F1,F2,F3])

    return F

def RichtMyer(q,dt,dx,Tfin,gamma):
    t = 0.0
    q = q
    while t < Tfin:
        qn = q.copy()
        # Create Flux Vector

        F = Flux_Vec(qn,gamma)

        # Start Predictor Step
        qp = np.roll(qn,-1)
        qm = np.roll(qn, 1)
        fp = np.roll(F,-1)
        fm = np.roll(F,1)

        qph = 0.5*(qp+q) - dt/(dx*2)*(fp-F)
        qmh = 0.5*(q+qm) - dt/(dx*2)*(F-fm)

        Fqph = Flux_Vec(qph,gamma)
        Fqmh = Flux_Vec(qmh,gamma)

        # Start Corrector Step
        q = qn - dt/dx*(Fqph - Fqmh)

        # Reinitialize Boundaries
        q[:,0] = qn[:,0]
        q[:,-1] = qn[:,-1]

        # Advance Time

        t = t+dt

    # PostProcess
    Density = q[0]
    Velocity = q[1]/Density
    E_fin = q[2]/Density
    Pressure = (gamma-1)*Density*(E_fin-0.5*Velocity**2)
    SpeedOfSound = np.sqrt(gamma*Pressure/Density)
    MachNumber = Velocity/SpeedOfSound

    return Density, Velocity, Pressure, SpeedOfSound, MachNumber