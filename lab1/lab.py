import matplotlib.pyplot as plt
import numpy as np

# All approximation-methods
def approximation1(x, h, func):
    return (func(x + h) - func(x)) / h
def approximation2(x, h, func):
    return (func(x) - func(x - h)) / h
def approximation3(x, h, func):
    return (func(x + h) - func(x - h)) / (2 * h)
def approximation4(x, h, func):
    return 4 / 3 * approximation3(x, h, func) - 1 / 3 * approximation3(x, 2 * h, func)
def approximation5(x, h, func):
    return 3 / 2 * approximation3(x, h, func) - 3 / 5 * approximation3(x, 2 * h, func) + \
         + 1 / 10 * approximation3(x, 3 * h, func)

# All functions for approximation
class func1:
    name = r'$\sin{x^{2}}$'
    @staticmethod
    def calc_func(x):
        return np.sin(x ** 2)
    @staticmethod
    def calc_deriv(x):
        return 2 * x * np.cos(x ** 2)

class func2:
    name = r'$\cos({\sin(x)})$'
    @staticmethod
    def calc_func(x):
        return np.cos(np.sin(x))
    @staticmethod
    def calc_deriv(x):
        return -np.sin(np.sin(x)) * np.cos(x)

class func3:
    name = r'$\exp({\sin({\cos(x)})})$'
    @staticmethod
    def calc_func(x):
        return np.exp((np.sin(np.cos(x))))
    @staticmethod
    def calc_deriv(x):
        return -np.sin(x) * np.cos(np.cos(x)) * np.exp(np.sin(np.cos(x)))

class func4:
    name = r'$\ln(x + 3)$'
    @staticmethod
    def calc_func(x):
        return np.log(x + 3)
    @staticmethod
    def calc_deriv(x):
        return 1 / (x + 3)

class func5:
    name = r'$(x + 3)^{0.5}$'
    @staticmethod
    def calc_func(x):
        return (x + 3) ** 0.5
    @staticmethod
    def calc_deriv(x):
        return 0.5 * (x + 3) ** (-0.5)


def calc_diff(x, h, approximation, func):
    return np.abs(approximation(x, h, func.calc_func) - func.calc_deriv(x))

steps = [(2 / 2**i) for i in range(1, 22)]
x = 10

approximations = [approximation1, approximation2, approximation3, approximation4, approximation5]
funcs = [func1, func2, func3, func4, func5]

for func, i in zip(funcs, range(1, 6)):
    y_list = []
    for method in approximations:
        y = [calc_diff(x, h, method, func) for h in steps]
        y_list.append(y)
    
    plt.title("function = " + func.name)
    plt.ylabel("Absolute error(log)")
    plt.xlabel("Step value(log)")
    
    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    plt.semilogx()
    plt.semilogy()
    
    plt.plot(steps, y_list[0], '-', steps, y_list[1], '--', steps, y_list[2], '-.', steps, \
             y_list[3], ':', steps, y_list[4])
    
    plt.savefig("graphs/funcGraph" + str(i) + ".png")
    plt.show()
