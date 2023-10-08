import numpy as np
import matplotlib.pyplot as plt
import copy 

class Matrix:
    def __init__(self, hight, width, elem_init):
        self.data = []
        self.hight = hight
        self.width = width

        for i in range(1, hight + 1):
            line = []
            for j in range(1, width + 1):
                line.append(elem_init(i, j))
            self.data.append(line)

    def dump(self):
        for line in self.data:
            for elem in line:
                print(format(elem, '.6f'), end=' ')
            print()
    
    def findMaxInLine(self, line_number, start_row = 0):
        line = self.data[line_number]
        idx_of_max = 0
        max_value = line[start_row];

        for i in range(start_row + 1, self.width):
            if line[i] > max_value and line[i] != 0:
                idx_of_max = i;
                max_value = line[i]

        return max_value, idx_of_max

    def findMaxInRow(self, row_number, start_line = 0):
        idx_of_max = start_line
        max_value = self.data[start_line][row_number];     

        for i in range(start_line + 1, self.hight):
            elem = self.data[i][row_number]
            if elem > max_value and elem != 0:
                idx_of_max = i;
                max_value = elem

        return max_value, idx_of_max

    # Method for supporting calculus of fixed accuracy
    # (round to 0 all elems less that 1e-6)
    def roundToZero(self):
        for i in range(self.hight):
            for j in range(self.width):
                if (abs(self.data[i][j] - 0) < 1e-6):
                    self.data[i][j] = 0

    def __mul__(self, column: list):
        res = []
        for i in range(self.hight):
            elem = 0
            for j in range(self.width):
                elem += self.data[i][j] * column[j]

            res.append(elem)

        return res

    #def transpose
    #def calc_det
    #def calc_own_values
    #def calc_own_vectors

# Gauss method
def solveGauss(matrix_in: Matrix, right_side_in: list):
    matrix = copy.deepcopy(matrix_in)
    right_side = copy.deepcopy(right_side_in)

    # Direct course
    for i in range(matrix.width):
        (max_row_elem, idx_of_max) = matrix.findMaxInRow(i, i)
        line_with_max_elem = matrix.data[idx_of_max]

        matrix.data[i], matrix.data[idx_of_max] = matrix.data[idx_of_max], matrix.data[i]
        right_side[i], right_side[idx_of_max] = right_side[idx_of_max], right_side[i]

        for j in range(i + 1, matrix.hight):
            curr_line = matrix.data[j];
            curr_line_elem_in_row = curr_line[i]

            coef = curr_line_elem_in_row / max_row_elem
            right_side[j] -= right_side[i] * coef
            
            for k in range(i, matrix.width):
                curr_line[k] -= line_with_max_elem[k] * coef

    matrix.roundToZero()
    
    # Reverse course 
    solutions = []
    for i in range(matrix.hight - 1, 0 - 1, -1):
        xi = right_side[i] / matrix.data[i][i]
        solutions.insert(0, xi)
        
        for j in range(i, 0 - 1, -1):
            right_side[j] -= matrix.data[j][i] * xi
            matrix.data[j][i] = 0
    
    return solutions
    
def solveJacobi(matrix_in: Matrix, right_side: list):
    # L- Lower triangular matrix
    L = copy.deepcopy(matrix_in)
    # D - Diagonal matrix
    D = copy.deepcopy(matrix_in)
    # U - Upper triangular matrix
    U = copy.deepcopy(matrix_in)
    
    for i in range(matrix_in.hight):
        for j in range(i, matrix_in.hight):
            U.data[i][j] = 0
            L.data[j][i] = 0

            if (i != j):
                D.data[i][j] = 0
                D.data[j][i] = 0

 
n = 3

matrix = Matrix(n, n, lambda i, j: 1 * (j == i) + (1 / (j + i)) * (i != j))
right_side = [1 / i for i in range(1, n + 1)]

resGauss = solveGauss(matrix, right_side)
resJacobi = solveJacobi(matrix, right_side)
