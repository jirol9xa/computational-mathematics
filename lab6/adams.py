import numpy as np
import random as rnd
import matplotlib.pyplot as plt

from numpy import array
from numpy import copy

def PrintHeader():
    plt.title("Phase trajectory")
    plt.ylabel("y")
    plt.xlabel("x")

def GenU(u0, rnd_size):
    return [u0[0] + rnd_size * rnd.randrange(-500, 500), u0[1] + rnd_size * rnd.randrange(-500, 500)]

def F(U):
    return [U[1], U[0] ** 2 - 1]


def NextU1(U, tau, f_arr):
    f_arr.append(F(U))
    return copy(U) + tau * (array(f_arr[-1]))


def NextU2(U, tau, f_arr):
    f_arr.append(F(U))

    last = array(f_arr[-1])
    penult = array(f_arr[-2])

    return copy(U) + tau * (1.5 * last - 0.5 * penult)


def NextU3(U, tau, f_arr):
    f_arr.append(F(U))

    last = array(f_arr[-1])
    penult = array(f_arr[-2])
    penult_prev = array(f_arr[-3])

    return copy(U) + tau * (23.0/12.0 * last - 16.0/12.0 * penult + 5.0/12.0 * penult_prev)



def AddRes(x_arr, y_arr, res):
    x_arr = np.append(x_arr, res[0])
    y_arr = np.append(y_arr, res[1])

def main():
    np.set_printoptions(floatmode='maxprec', suppress=True)

    T = 1e-3
    rnd_size = 0.001

    plt_rad = 100
    u0 = [1.0, 0.0]
    points_amnt = 800

    for point in range(points_amnt):
        f_arr = []
        u = GenU(u0, rnd_size)
        x_arr = array([u[0]])
        y_arr = array([u[1]])

        res = u
        for func in [NextU1, NextU2]:
            last = copy(res)
            res = func(last, T, f_arr)
            AddRes(x_arr, y_arr, res) 

        for t in range(600):
            last = np.copy(res)
            res = NextU3(last, T, f_arr)
            if (res[0] * res[0] + res[1] * res[1]) > plt_rad:
                break
            AddRes(x_arr, y_arr, res)

        print(f"\rPlotting point = {point}", end="")
        plt.plot(x_arr, y_arr, '.-', ms=1)
    plt.grid()

    print("\nDone")

    PrintHeader()
    plt.show()

if __name__ == "__main__":
    main()
