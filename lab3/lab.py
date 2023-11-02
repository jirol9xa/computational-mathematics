import math
import matplotlib.pylab as plt

eps = 1e-6

def newton_method(x, func, deriv_func, solutions_chain: list = []):
    solutions_chain.append(x)
    if func(x) < eps:
        return solutions_chain
    next_x = x - func(x) / deriv_func(x)
    return newton_method(next_x, func, deriv_func, solutions_chain)


def simple_iteration_method(x, func, next_x_func, solutions_chain: list = []):
    solutions_chain.append(x)
    if func(x) < eps:
        return solutions_chain
    next_x = next_x_func(x)
    return simple_iteration_method(next_x, func, next_x_func, solutions_chain)


def func(x):
    print("x = ", x)
    return 2 * x**2 + 5 * x - 3


def deriv_func(x):
    return 4 * x + 5


def next_x_func(x):
    return  (3 - 2 * x**2) / 5 #math.sqrt((3 - 5 * x) / 2)


def find_nearest_solution_newton(x):
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



find_nearest_solution_newton(4)
find_nearest_solution_newton(-5)

# Can't solve prev equasion with simple iter method, need to pick another one
find_nearest_solution_iterations()
find_nearest_solution_iterations()
