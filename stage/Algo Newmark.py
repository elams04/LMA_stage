import numpy as np
from scipy.sparse import csr_matrix

def Newmark(obj, F, m, c, k, t, implementation):
    # Newmark
    # Solver for the EDP integration

    if c is None or c.size == 0:
        damping = 0
        c = csr_matrix(np.zeros(k.shape))

    # time def
    dt = t[1] - t[0]  # size of the time step
    nt = int((t[-1] - t[0]) // dt)  # number of time steps
    sdof = k.shape[0]

    # initializing the displacement, velocity and acceleration matrices
    depl = np.zeros((sdof, nt))
    vel = np.zeros((sdof, nt))
    accl = np.zeros((sdof, nt))
    Reff = np.zeros((sdof, nt))

    # Parameters for Newmark time integration
    beta = obj.param[0]
    gamma = obj.param[1]

    if implementation == 'Bathe':
        if damping == 0:
            # Solve for initial accelerations
            accl[:, 0] = np.linalg.solve(m, F[:, 0] - k @ depl[:, 0])
            # Calculating integration constants
            a0 = 1 / (beta * dt**2)
            a1 = gamma / (beta * dt)
            a2 = 1 / (beta * dt)
            a3 = (1 / (2 * beta)) - 1
            a4 = (gamma / beta) - 1
            a5 = (dt / 2) * (gamma / beta - 2)
            a6 = dt * (1 - gamma)
            a7 = gamma * dt
            keff = k + a0 * m
            # time step starts
            for it in range(nt):
                Reff[:, it] = F[:, it + 1] + m @ (a0 * depl[:, it] + a2 * vel[:, it] + a3 * accl[:, it])
                # solving for displacements at time (it+dt)
                depl[:, it + 1] = np.linalg.solve(keff, Reff[:, it])
                # calculating velocities and accelerations at time (it+dt)
                accl[:, it + 1] = a0 * (depl[:, it + 1] - depl[:, it]) - a2 * vel[:, it] - a3 * accl[:, it]
                vel[:, it + 1] = vel[:, it] + a6 * accl[:, it] + a7 * accl[:, it + 1]
        else:
            # initial accelerations
            accl[:, 0] = np.linalg.solve(m, F[:, 0] - c @ vel[:, 0] - k @ depl[:, 0])
            # constants
            a0 = 1 / (beta * dt**2)
            a1 = gamma / (beta * dt)
            a2 = 1 / (beta * dt)
            a3 = (1 / (2 * beta)) - 1
            a4 = (gamma / beta) - 1
            a5 = (dt / 2) * (gamma / beta - 2)
            a6 = dt * (1 - gamma)
            a7 = gamma * dt

            keff = k + a0 * m + a1 * c
            for it in range(nt):
                Reff[:, it] = F[:, it + 1] + m @ (a0 * depl[:, it] + a2 * vel[:, it] + a3 * accl[:, it]) + c @ (a1 * depl[:, it] + a4 * vel[:, it] + a5 * accl[:, it])
                # solving displacements at time (it+dt)
                depl[:, it + 1] = np.linalg.solve(keff, Reff[:, it])
                # velocities and accelerations at time (it+dt)
                accl[:, it + 1] = a0 * (depl[:, it + 1] - depl[:, it]) - a2 * vel[:, it] - a3 * accl[:, it]
                vel[:, it + 1] = vel[:, it] + a6 * accl[:, it] + a7 * accl[:, it + 1]
    elif implementation == 'a-Form-Hughes':
        raise NotImplementedError("a-Form-Hughes implementation is not provided.")

