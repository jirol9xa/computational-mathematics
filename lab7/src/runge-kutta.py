import numpy as np
from numpy import array, zeros

class RungeKutta:
  def __init__(self, matr, t0, y0, f):
    self.size_ = len(matr) - 1
    self.a_, self.B, self.C = self.GetCoeffs(matr) 
    self.explicit_ = self.isExplicit(matr)
    self.y_ = y0
    self.dy_ = f(t0, y0)
    self.time_ = t0
    self.f_ = f

    self.y_trace = []
    self.dy_trace = []
    for i in range(0, len(self.y_)): 
      self.y_trace.append([]) 
      self.dy_trace.append([])
    self.time_trace = []

    self.SaveTrace()


  def GetCoeffs(self, matr):
    return [matr[:self.size_, 1:],               \
            matr[self.size_, 1:].tolist()[0],    \
            matr[:self.size_, 0].T.tolist()[0]]


  def isExplicit(self, matr):
    return all(matr[i][j] == 0 for i in range(len(matr)) for j in range(i + 1, len(matr[i])))

  def SaveTrace(self):
    for i in range(0, len(self.y_)):
      (self.y_trace[i]).append(self.y_[i])
      (self.dy_trace[i]).append(self.dy_[i])

    self.time_trace.append(self.time_)

  def NextIteration(self, dt, explicit=None):
    if dt < 0:
      raise RuntimeError("Time can not be negative value")

    if explicit:
      self.ExplicitNextIteration(dt)
    else:
      self.ImplicitNextStep(dt)

    self.time_ += dt
    self.dy_ = self.f_(self.time_, self.y_)
    self.SaveTrace()


  def GetIterAmnt(self, step, finish_time):
    return len(np.arange(self.time_, finish_time + step, step))
  

  def ExplicitNextIteration(self, dt):
    if dt < 0:
      raise RuntimeError("Time can not be negative value")

    # Get K for every single iter
    K = np.array(np.zeros(self.size_))
    for i in range(0, self.size_):
      t_s = self.C[i] * dt
      y_arr = np.array(np.zeros(len(self.y_)))
      y_arr = self.y_
      for j in range(0, i):
        y_arr += self.a_[i, j] * dt * K[j]
      K[i] = self.f_(t_s, y_arr)
    
    # Calc y(n + 1)
    for i in range(0, self.size_):
      self.y_ += dt * self.B[i] * np.copy(K[i]) 


  def ImplicitNextStep(self, dt):
    if dt < 0:
      raise RuntimeError("Time can not be negative value")

    
    Equation = lambda y_np1: y_np1 - dt * self.f_(self.time_, y_np1) - self.y_
    self.y_ = Solve(Equation, x0=np.copy(self.y_))


  def Start(self, step, finish_time, explicit=None):
      for _ in range(1, self.GetIterAmnt(step, finish_time)):
        self.NextIteration(step, explicit=explicit)
