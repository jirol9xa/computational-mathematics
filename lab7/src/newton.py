import numpy as np


def DiffV5(func, h, x, xid=0):
  h_arr = np.zeros(len(x))
  h_arr[xid] = h
  return 3 * (func(x + h_arr) - func(x - h_arr)) / (2 * 2 * h) - \
         3 * (func(x + 2 * h_arr) - func(x - 2 * h_arr)) / (5 * 4 * h) + \
         1 * (func(x + 3 * h_arr) - func(x - 3 * h_arr)) / (10 * 6 * h)


def CalcJacobian(F, x):
  h = 10**(-5)
  SizeX = len(x)
  res = np.zeros([SizeX, SizeX])
  for fid in range(0, SizeX):
    for xid in range(0, SizeX):
      func = lambda x: F(x)[fid]
      res[fid, xid] = DiffV5(func, h, x, xid)

  return res


def Solve(func, x0, eps=10**-4):
  sol = x0
  cnt = 0
  J = CalcJacobian(func, sol)
  J = np.linalg.inv(J)

  while (cnt <= 9):
    F = np.matrix(func(sol))
    sol = np.array((np.matrix(sol).T - J * F.T).T)[0]
    cnt += 1

  return sol