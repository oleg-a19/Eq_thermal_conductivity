import numpy as np



class thermal_conductivity2D:
    def __init__(self, grid, T):
        self.grid = grid
        self.T = T
        self.u = np.empty(len(grid.points))

        self.cp = np.empty(len(grid.points))
        self.cp.fill(840.)
        self.rho = np.empty(len(grid.points))
        self.rho.fill(37.)
        self.lamb = np.empty(len(grid.points))
        self.lamb.fill(0.033)

    def set_init_condit (self, T_init):
        self.u.fill(T_init)

    def set_boundary_condit(self, uL, uR):
        for k in range(self.grid.K):
            for m in range(self.grid.M):
                if k == 0:
                    i = m+k*self.grid.M
                    self.u[i] = 1./(uL[2]-uL[0]/uL[1]) * (uL[3]-uL[0]/uL[1]*self.u[m+(k+1)*self.grid.M])
                elif k == self.grid.K-1:
                    i = m+k*self.grid.M
                    self.u[i] = 1./(uR[2]-uR[0]/uR[1]) * (uR[3]-uR[0]/uR[1]*self.u[m+(k-1)*self.grid.M])
                elif m == self.grid.M-1:
                    i = m+k*self.grid.M
                    self.u[i] = self.u[i-1]
                elif m == 0:
                    i = m+k*self.grid.M
                    self.u[i] = self.u[i+1]

    def solve(self):
        u0 = np.copy(self.u)
        hx = (self.grid.corner_points[3].x-self.grid.corner_points[0].x)/self.grid.K
        hy = (self.grid.corner_points[1].y-self.grid.corner_points[0].y)/self.grid.M
        h = np.array([hx, hy])
        self.set_boundary_condit((-0.033,h[0],0., 4e2), (0.033,h[0],0., 0.))
        t = 0.
        while t < self.T:
            dt = 0.24*np.min(self.cp)*np.min(self.rho)*np.min(h)**2 / np.max(self.lamb)
            for k in range(1, self.grid.K-1):
                for m in range(1, self.grid.M-1):
                    # mk, m+1k, mk+1, m-1k, mk-1, m+1/2k, mk+1/2, m-1/2k, mk-1/2
                    # 0   1     2     3     4     0и1     0и2     0и3     0и4     
                    i = np.array([m+k*self.grid.M, (m+1)+k*self.grid.M, m+(k+1)*self.grid.M,
                        (m-1)+k*self.grid.M, m+(k-1)*self.grid.M])
                    sig = np.array([dt / h[0] / self.rho[i[0]] / self.cp[i[0]], dt / h[1] / self.rho[i[0]] / self.cp[i[0]]])

                    self.u[i[0]] = u0[i[0]] + (sig[0]*((self.lamb[i[0]]+self.lamb[i[1]])*0.5*(u0[i[1]]-u0[i[0]])/h[0]
                                                    -(self.lamb[i[0]]+self.lamb[i[3]])*0.5*(u0[i[0]]-u0[i[3]])/h[0])
                                            +
                                            sig[1]*((self.lamb[i[0]]+self.lamb[i[2]])*0.5*(u0[i[2]]-u0[i[0]])/h[1]
                                                    -(self.lamb[i[0]]+self.lamb[i[4]])*0.5*(u0[i[0]]-u0[i[4]])/h[1]))
            self.set_boundary_condit((-0.033,h[0],0, 4e2), (0.033,h[0],0, 0.))
            u0 = np.copy(self.u)
            t += dt
            if t/dt%1000 - 0 <= 1e-3 :
                print(dt, t)