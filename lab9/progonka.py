import numpy as np
import matplotlib.pyplot as plt

h = 0.005
N = int(1 / h)

def f(x):
    return np.cos(2*np.pi*x)

def P_squared(x):
    return 10 + np.sin(2*np.pi*x)

def a(k):
    return 1

def b(k):
    return 2 + P_squared(h*k) * h**2
    
def c(k):
    return 1

def phi(k):
    return f(h*k) * h**2

coeffs1 = []
alpha = c(0) / b(0)
beta = - phi(0) / b(0)
gama =  a(0) / b(0)

coeffs1 = np.zeros(shape=(N + 1, 3))
coeffs1[1] = [alpha, beta, gama]

print(alpha, beta, gama)

for i in range (1, N):
    alpha = c(i) / (b(i) - coeffs1[i][0]*a(i))
    beta  = (a(i)*coeffs1[i][1] - phi(i)) / (b(i) - coeffs1[i][0]*a(i))
    gama  = a(i)*coeffs1[i][2] / (b(i) - coeffs1[i][0]*a(i))
    coeffs1[i + 1] = [alpha, beta, gama]

# Second iteration
mu = np.zeros(N + 1)
nu = np.zeros(N + 1)

mu[N] = -c(N) / (a(N)*(coeffs1[N][0] + coeffs1[N][2]) - b(N))
nu[N] = (phi(N) - a(N)*coeffs1[N][1]) / (a(N)*(coeffs1[N][0] + coeffs1[N][2]) - b(N))

for i in range (N , 0, -1):
    mu[i - 1] = coeffs1[i][0] * mu[i] + coeffs1[i][2] * mu[N]
    nu[i - 1] = coeffs1[i][1] + coeffs1[i][0] * nu[i] + coeffs1[i][2]*nu[N]


y = np.zeros(N + 1)

y_0 = nu[0] / (1 - mu[0])
y[N] = mu[N] * y_0 + nu[N]

#Third iteration and get solution
for i in range(N, 0, -1):
    y[i - 1] = coeffs1[i][0]*(mu[i]*y_0 + nu[i]) + coeffs1[i][1] + coeffs1[i][2]*(mu[N]*y_0 + nu[N])

x = np.arange(0, 1 + h, h)

plt.title("Solution")
plt.grid()
plt.plot(x, y)
plt.show()
