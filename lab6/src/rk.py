import matplotlib.pyplot as plt
import numpy as np
import random as rnd

from numpy import array, copy

def NextU(U, T, b, k_arr):
    res = [0, 0]

    for i in range(len(b)):
        res += k_arr[i] * b[i] * T

    return copy(U) + res


def PrintHeader():
    plt.title("Phase trajectory")
    plt.ylabel("y")
    plt.xlabel("x")


def F(U):
    return [U[1], U[0] * U[0] - 1]


def K(U, s, tau, A):
    if s <= 1:
        return F(U)

    k_arr = np.empty((s + 1, 2))
    k_arr[0] = F(U)

    for i in range(1, s + 1):
        U_new = copy(U)
        for j in range(i):
            U_new += tau * A[i - 1][j] * k_arr[j]
        k_arr[i] = F(U_new)

    return k_arr


def main():
    np.set_printoptions(floatmode='maxprec', suppress=True)

    plot_radius = 100
    u0 = [-1.0, 0.0]
    points_amnt = 800

    A = np.empty((3, 3))
    A = [[0.5, 0,   0],
         [0,   0,5, 0],
         [0,   0,   1]]

    T = 0.01
    rnd_size = 0.001

    b = array([1.0/6, 1.0/3, 1.0/3, 1.0/6])

    for point in range(points_amnt):
        u = [u0[0] + rnd_size * rnd.randrange(-5, 5000), u0[1] + rnd_size * rnd.randrange(-200, 200)]
        x_arr = array([u[0]])
        y_arr = array([u[1]])

        for t in range(0, 100):
            last_u = [x_arr[-1], y_arr[-1]]
            k_arr = K(last_u, 3, T, A)
            res = NextU(last_u, T, b, k_arr)
            if (res[0] * res[0] + res[1] * res[1]) > plot_radius:
                break
            x_arr = np.append(x_arr, res[0])
            y_arr = np.append(y_arr, res[1])

        plt.plot(x_arr, y_arr, '--', ms=1)

    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    plt.grid()

    print("\nDone")

    plt.show()

if __name__ == "__main__":
    main()
