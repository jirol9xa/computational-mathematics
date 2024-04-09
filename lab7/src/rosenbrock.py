import numpy as np

class RosenBrok():
    def __init__(self, f, t, x, h):
        self.x_, self.f_, self.t_, self.h_ = x, f, t, h
    
    def RosenBrokMethod(self, f, x_0, t_1, t_2, N):
        dim = len(x_0)
        x = np.empty(shape=(N + 1, dim))
        x[0] = x_0

        for n in range(N):
            t = np.linspace(t_1, t_2, num = N + 1, endpoint = True)
            dt = (t_2 - t_1) / N
            
            J = self.JacMatr(f, t[n], x[n], dt)
            M = np.eye(dim) - (0.5 + 0.5 * 1j) * dt * J
            f = np.empty(dim)
            for i in range(dim):
                f[i] = f[i](t[n] + dt * 0.5, x[n])

            w = np.linalg.solve(M, f)
            x[n + 1] = x[n] + dt * w.real

        return t, x

    def JacMatr(self, f, t, x):
        J = np.empty(shape=(len(f), len(x)))
        for i in range(len(f)):
            for j in range(len(x)):
                J[i, j] = self.PartDeriv(self.f_[i],t, self.x_, self.h_, j)

        return J
    
    def PartDeriv(self, j):
        H = np.zeros_like(x)
        H[j] = self.h_
        return 3 * (self.f_(self.t_, self.x_ + H) - self.f_(self.t_, self.x_ - H)) / (4 * self.h_) -       \
               3 * (self.f_(self.t_, self.x_ + 2*H) - self.f_(self.t_, self.x_ - 2*H)) / (20 * self.h_) +  \
               (self.f_(self.t_, self.x_ + 3*H) - self.f_(self.t_, self.x_ - 3*H)) / (60 * self.h_)

