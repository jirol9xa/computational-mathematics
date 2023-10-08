import numpy as np
import matplotlib.pyplot as plt

class Matrix:
    data = []
    hight = 0
    width = 0

    def __init__(self, hight, width, elem_init):
        self.hight = hight
        self.width = width

        for i in range(1, hight + 1):
            line = []
            for j in range(1, width + 1):
                line.append(elem_init(i, j))
            self.data.append(line)

    def dump(self):
        for line in self.data:
            print(line)
    
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

    #def transpose
    #def calc_det
    #def calc_own_values
    #def calc_own_vectors


# Gauss method
def solveGauss(matrix: Matrix):
    for i in range(matrix.width):
        (max_row_elem, idx_of_max) = matrix.findMaxInRow(i, i)
        line_with_max_elem = matrix.data[idx_of_max]

        for j in range(i + 1, matrix.hight):
            curr_line = matrix.data[j];
            curr_line_elem_in_row = curr_line[i]

            for k in range(i, matrix.width):
                curr_line[k] -= line_with_max_elem[k] * (curr_line_elem_in_row / max_row_elem)

n = 3
matrix = Matrix(n, n, lambda i, j: 1 * (j == i) + (1 / (j + i)) * (i != j))
solveGauss(matrix)
matrix.roundToZero()
