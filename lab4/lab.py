import matplotlib.pylab as plt
import numpy as np


data = {    
            1910 : 92228496, 
            1920 : 106021537,
            1930 :123202624, 
            1940 : 132164569,
            1950 : 151325798, 
            1960 : 179323175, 
            1970 : 203211926,
            1980 : 226545805,
            1990 : 248709873, 
            2000 : 281421906
        }

amount = list(data.values())
years = list(data.keys())


class NewtonAprox:
    def __init__(self, years, amount):
        self.b_coef = [[0 for _ in range(len(years))] for _ in range(len(years))]
        self.years = years
        self.amount = amount
        self.name = "Newton"

        for i in range(1, len(years)):
            for j in range(len(years) - i):
                if (i == 1):
                    self.b_coef[i][j] = (amount[j + 1] - amount[j]) / (years[j + 1] - years[j])
                else:
                    self.b_coef[i][j] = (self.b_coef[i - 1][j + 1] - self.b_coef[i - 1][j]) / (years[j + i] - years[j])

    def getName(self):
        return self.name

    def calcValue(self, x):
        res = self.amount[0]

        for i in range(1, len(years)):
            tmp_res = self.b_coef[i][0] 
            for j in range(1, i + 1):
                tmp_res *= (x - self.years[j - 1])
            res += tmp_res

        print("x = ", x, "value = ", res)
        return res


class SplineAprox:
    def __init__(self, years, amount):
        self.years = years
        self.amount = amount
        self.name = "Spline"    

    def getName(self):
        return self.name

    def calcValue(self, x):


methods = [NewtonAprox(years, amount), SplineAprox(years, amount)]
years_plot = [i for i in range(1910, 2020, 10)]

for i in range(len(methods)):
    population = [methods[i].calcValue(j) for j in years_plot]

    plt.title(methods[i].getName() + " method for population 1910-2010")
    plt.ylabel("Population")
    plt.xlabel("Year")

    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')

    plt.plot(years_plot, population)
    plt.savefig("images/" + methods[i].getName())
    plt.show()

    print("f(2010) = ", methods[i].calcValue(2010))
