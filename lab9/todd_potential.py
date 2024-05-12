import numpy as np
import matplotlib.pyplot as plt
import sys

delta_tau = 0.001
t_max = 120
A_delta = 0.01
f_delta = 0.001

def f(y_n):
    return np.array([y_n[1], 1 - np.exp(y_n[0])])

def k_1(y_n):
    return f(y_n)

def k_2(y_n, k_1):
    y = y_n + k_1 * delta_tau / 2
    return f(y)

def k_3(y_n, k_2):
    y = y_n + k_2 * delta_tau / 2
    return f(y)

def k_4(y_n, k_3):
    y = y_n + k_3 * delta_tau
    return f(y)

def solution_periods(T : int) -> list:
    periods = set()

    for i in range (1, int(T**0.5) + 1):
        if (T % i == 0):
            periods.add(i)
            periods.add(T // i)

    return list(periods)

periods = solution_periods(t_max)


def plot_solution(x : np.ndarray):
    plt.xlabel("t, time")
    plt.ylabel("x, coordinate")
    plt.title("Solution")
    plt.grid()

    t = np.arange(0, t_max, delta_tau)

    plt.plot(t, x)

def right_correct(x, A) -> bool:
    deriv_left  = A - A_delta
    deriv_right = A + A_delta

    deriv = (x[len(x) - 1] - x[len(x) - 2]) / delta_tau

    if (deriv > deriv_left and deriv < deriv_right):
        return True
    
    return False


# to find f(alpha) we will use Runge-Kutte method
def solve_koshi(x_start, v_start, t_end):
    y = []
    y_0 = np.array([x_start, v_start])
    y.append(y_0)

    for i in range (1, int(t_end / delta_tau)):
        k_1_ = k_1(y[i - 1])
        k_2_ = k_2(y[i - 1], k_1_)
        k_3_ = k_3(y[i - 1], k_2_)
        k_4_ = k_4(y[i - 1], k_3_)

        y_i = y[i - 1] + delta_tau*(k_1_ + 2*k_2_ + 2*k_3_ + k_4_) / 6
        y.append(y_i)

    #separate solutions:
    x = []
    v = []
    for i in range (0, len(y)):
        x.append(y[i][0])
        v.append(y[i][1])

    return x, v

def solve_fire(x_start, x_end, t_end, A):
    alpha_left  = A - A_delta
    alpha_right = A + A_delta

    h = A_delta / 1000

    alpha, x, v = find_root(x_start, x_end, t_end, alpha_left, alpha_right, h)

    if not ((alpha > alpha_left) and (alpha < alpha_right)):
        return -1, -1

    return x, v

def find_root(x_start, x_end, t_end, left, right, h):
    alpha = (left + right) / 2

    x = []
    v = []
    x, v = solve_koshi(x_start, alpha, t_end)
    f_alpha = x[len(x) - 1] - x_end
    x, v = solve_koshi(x_start, alpha + h, t_end)
    f_alpha_h = x[len(x) - 1] - x_end

    while (abs(f_alpha) > f_delta):
        f_alpha_der = (f_alpha_h - f_alpha) / h
        alpha = alpha - f_alpha / f_alpha_der

        x, v = solve_koshi(x_start, alpha, t_end)
        f_alpha = x[len(x) - 1] - x_end
        x, v = solve_koshi(x_start, alpha + h, t_end)
        f_alpha_h = x[len(x) - 1] - x_end

    return alpha, x, v

def extend_periodically(x):
    T = int(len(x)*delta_tau)

    num_periods = t_max // T
    x_extended = []

    for period in range(0, num_periods):
        for elem in x:
            x_extended.append(elem)

    return np.array(x_extended)


if __name__ == "__main__":
    solutions = []
    A = np.arange(10, 100, 0.1)

    a = 14.52
    x_all = []
    for T in periods:
        x, v = solve_fire(0, 0, T, a)
        if x != -1 and right_correct(x, a):
            print("Solution Found!!!, a =", a)
            x_all.append(extend_periodically(x))
        solutions.append([a, x_all])
