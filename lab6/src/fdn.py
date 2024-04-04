import matplotlib.pyplot as plt
import numpy as np
import random as rnd

from numpy import array
from numpy import copy

def BackDIff(u_arr, p):
    if p == 1:
        return u_arr[-1] - u_arr[-2]
    return BackDIff(u_arr, p - 1) - BackDIff(u_arr[:-1], p - 1)


def NumDiffArr(func, x, h, step):
    return (func(x + 3 * h) - func(x - 3 * h)) / 6.0 / step


def J(U, func):
    size = len(U)
    res = np.zeros((size, size))
    step = 1e-5
    h = np.zeros(size)

    for i in range(size):
        h[i] = step
        res[i] = NumDiffArr(func, U, h, step)
        h[i] = 0

    return res


def F(U):
    return array([U[1], U[0] * U[0] - 1])


def NextU1(U, tau, unused):
    return copy(U) + tau * (array(F(U)))


def NextU2(U, tau, u_arr):
    return copy(U) + tau * (3.0/2.0 * F(u_arr[-1]) - 1.0/2.0 * F(u_arr[-2]))


def NextU3(U, tau, u_arr):
    return copy(U) + tau * (23.0/12.0 * F(U) - 16.0/12.0 * F(u_arr[-1]) + 5.0/12.0 * F(u_arr[-2]))


def FDNNextU(alpha_arr, beta, tau, func, u_arr):
    last = u_arr[-1]
    delta = (alpha_arr[2] * last + alpha_arr[1] * u_arr[-2] + alpha_arr[0] * u_arr[-3])
    gamma = tau * beta
    new = gamma * F(last) + delta

    while (np.linalg.norm(last - new) > 1e-4):
        last = new
        new = gamma * F(last) + delta

    return new


def PrintHeader():
    plt.title("Phase trajectory")
    plt.ylabel("y")
    plt.xlabel("x")


def main():
    np.set_printoptions(floatmode='maxprec', suppress=True)
    tau = 1e-3
    rnd_size = 0.01

    plt_rad = 100
    u0 = [1.0, 0.0]

    alpha_arr = [2.0/11.0, -9.0/11.0, 18.0/11.0]
    beta = 6.0/11.0

    points_amnt = 600

    for _ in range(points_amnt):
        u = [u0[0] + rnd_size * rnd.randrange(-5, 500), u0[1] + rnd_size * rnd.randrange(-200, 200)]
        u_arr = [u]

        res = u
        for func in [NextU1, NextU2, NextU3]:
            last_u = copy(res)
            res = func(last_u, tau, u_arr)
            u_arr.append(res)   

        for _ in range(900):
            res = FDNNextU(alpha_arr, beta, tau, F, u_arr)
            if (res[0] * res[0] + res[1] * res[1]) > plt_rad:
                break
            u_arr.append(res)

        x_arr, y_arr = zip(*u_arr)
        plt.plot(x_arr, y_arr, '.-', ms=1)
    plt.grid()

    print("\nDone")

    plt.show()

if __name__ == "__main__":
    main()
