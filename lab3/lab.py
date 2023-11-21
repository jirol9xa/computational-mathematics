import math
import matplotlib.pylab as plt
import numpy as np

eps = 1e-6

def newton_method(x, func, deriv_func, solutions_chain: list = []):
    solutions_chain.append(x)
    if abs(func(x)) < eps:
        return solutions_chain
    next_x = x - func(x) / deriv_func(x)
    return newton_method(next_x, func, deriv_func, solutions_chain)


def simple_iteration_method(x, func, next_x_func, solutions_chain: list = []):
    solutions_chain.append(x)
    if abs(func(x)) < eps:
        return solutions_chain
    next_x = next_x_func(x)
    return simple_iteration_method(next_x, func, next_x_func, solutions_chain)


def newton_method_systems(x, func, jacobi_matr, solutions_chain: list = []):
    solutions_chain.append(x[0])
    if (abs(func(x)[0] ** 2 + func(x)[1] ** 2) < eps):
        return solutions_chain
    next_x = x - np.matmul(jacobi_matr(x), func(x))
    return newton_method_systems(next_x, func, jacobi_matr, solutions_chain)


def simple_iteration_method_systems(x, func, next_x_func, solutions_chain):
    solutions_chain.append(x[0])
    if (abs(func(x)[0] ** 2 + func(x)[1] ** 2) < eps):
        return solutions_chain
    next_x = next_x_func(x)
    return simple_iteration_method_systems(next_x, func, next_x_func, solutions_chain)


def find_nearest_solution_newton(x):
    def func(x):
        return 2 * x**2 + 5 * x - 3

    def deriv_func(x):
        return 4 * x + 5

    solutions = newton_method(x, func, deriv_func, [])
    iterations = [i + 1 for i in range(len(solutions))]

    plt.title("Newton method for x_0 = " + str(x))
    plt.ylabel("Solution value")
    plt.xlabel("Iterations")

    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')

    plt.plot(iterations, solutions)
    plt.savefig("images/Newton_" + str(x) + ".png")
    plt.show()

def find_nearest_solution_newton_system(x):
    jacobi_inv = lambda x: 1/(-np.sin(x[0] - 1) * np.sin(x[1]) -1) * np.array([[np.sin(x[1]), -1], [-1, -np.sin(x[0] - 1)]])
    func = lambda x: [ np.cos(x[0] - 1) + x[1] - 0.5, x[0] - np.cos(x[1]) - 3]

    solutions = newton_method_systems(x, func, jacobi_inv, [])
    iterations = [i for i in range(len(solutions))]

    print(len(solutions))

    plt.title("Newton method for system with (x_0, y_0) = " + str(x))
    plt.ylabel("Solution value")
    plt.xlabel("Iterations")

    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')

    plt.plot(iterations, solutions)
    plt.savefig("images/Newton_" + str(x) + ".png")
    plt.show()

def find_nearest_solution_iterations(x):
    def func(x):
        return x * 2**x - 1
    
    def next_x_func(x):
        return 1 / 2**x

    solutions = simple_iteration_method(x, func, next_x_func, [])
    iterations = [i + 1 for i in range(len(solutions))]

    plt.title("Simple iterations method for x_0 = " + str(x))
    plt.ylabel("Solution value")
    plt.xlabel("Iterations")

    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')

    plt.plot(iterations, solutions)
    plt.savefig("images/Simple_iterations_" + str(x) + ".png")
    plt.show()

def find_nearest_solution_iterations_systems(x):
    func = lambda x: [np.cos(x[0] - 1) + x[1] - 0.5, x[0] - np.cos(x[1]) - 3]
    next_x_func = lambda x: [3 + np.cos(x[1]), 0.5 - np.cos(x[0] - 1)]

    solutions = simple_iteration_method_systems(x, func, next_x_func, [])
    iterations = [i + 1 for i in range(len(solutions))]

    plt.title("Simple iterations method for x_0 = " + str(x))
    plt.ylabel("Solution value")
    plt.xlabel("Iterations")

    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')

    plt.plot(iterations, solutions)
    plt.savefig("images/Simple_iterations_" + str(x) + ".png")
    plt.show()



find_nearest_solution_newton(4)
find_nearest_solution_newton(-5)

find_nearest_solution_iterations(-0.1)

find_nearest_solution_newton_system([5, 0])
find_nearest_solution_iterations_systems([5, 0])
