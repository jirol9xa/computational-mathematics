import math
import matplotlib.pyplot as plt
import numpy as np

from numpy.linalg import inv
from math import sin,cos


def PrintHeader(method_name):
    plt.title("X vs Y. " + method_name)
    plt.xlabel("X")
    plt.ylabel("Y")


def MPIX(x, y):
    return 0.5 - cos(y - 2)


def MPIY(x, y):
    return sin(x + 2) - 1.5


def J(x: float, y:float):
    ret = np.empty((2, 2))
    ret = [[cos(x + 2), -1],
           [1,          -sin(y - 2)]]
    return ret


def F(x: float, y: float):
    return [sin(x + 2) - 1.5 - y, cos(y - 2) + x - 0.5]


def CalcNext(x_last, y_last, dF, F):
    tmp = (inv(dF) @ F(x_last, y_last))
    return [x_last - tmp[0], y_last - tmp[1]]


def main():
    eps = 0.001

    # MPI
    print("MPI:")

    x_last = 10.0
    y_last = 10.0

    x_arr = [x_last]
    y_arr = [y_last]

    x_new = MPIX(x_last, y_last)
    y_new = MPIY(x_last, y_last)
    while (abs(x_new - x_last) > eps) and (abs(y_new - y_last) > eps):
        print(x_new, y_new)
        x_last, y_last = x_new, y_new

        x_arr.append(x_last)
        y_arr.append(y_last)

        x_new, y_new = MPIX(x_last, y_last), MPIY(x_last, y_last)

    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    plt.grid()

    PrintHeader("MPI")

    plt.plot(x_arr, y_arr, '.-', ms=8.0)
    plt.savefig("img/system_mpi.png")
    plt.clf()

    #Newton
    print("Newton:")
    x_last = 5.0
    y_last = 5.0

    x_arr = [x_last]
    y_arr = [y_last]

    print(x_last, y_last)
    dF = J(x_last, y_last)
    x_new, y_new = CalcNext(x_last, y_last, dF, F)
    
    x_arr.append(x_new)
    y_arr.append(y_new)

    while (abs(x_new - x_last) > eps) and (abs(y_new - y_last) > eps):
        print(x_new, y_new)
        x_last, y_last = x_new, y_new

        dF = J(x_last, y_last)
        print(dF)
        print(inv(dF))
        x_new, y_new = CalcNext(x_last, y_last, dF, F)
        
        x_arr.append(x_new)
        y_arr.append(y_new)
        
    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    plt.grid()

    PrintHeader("Newthon")

    plt.plot(x_arr, y_arr, '--', ms=5.0)
    plt.savefig("img/system_newton.png")
    plt.clf()

if __name__ == '__main__':
    main()
