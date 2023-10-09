import numpy as np
import matplotlib.pyplot as plt
import copy 

class Matrix:
    def __init__(self, *args):
        if len(args) > 1:
            self.data = []
            self.hight = args[0]
            self.width = args[1]

            for i in range(1, self.hight + 1):
                line = []
                for j in range(1, self.width + 1):
                    line.append(args[2](i, j))
                self.data.append(line)

        elif len(args) == 1:
            self.data = args[0]
            self.hight = len(self.data)
            self.width = len(self.data[0])
        else:
            raise ValueError("Can't init matrix with that pack or args")


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


    def __add__(self, other):
        if self.hight != other.hight or self.width != other.width:
            raise ValueError("It is not possible to add matrices with different sizes")
        
        res = copy.deepcopy(self)
        for i in range(self.hight):
            for j in range(self.width):
                res.data[i][j] += other.data[i][j]

        return res


    def __mul__(self, other):
        if isinstance(other, list):
            res = []
            for i in range(self.hight):
                elem = 0
                for j in range(self.width):
                    elem += self.data[i][j] * other[j]

                res.append(elem)

            return res
        
        elif isinstance(other, int):
            res = copy.deepcopy(self);
            for i in range(res.hight):
                for j in range(res.width):
                    res.data[i][j] *= other

            return res

        if self.width != other.hight:
            raise ValueError("Can calculate the product of matrices only if \
                             the number of columns of the first matrix is equal \
                             to the number of rows of the second")
        
        res = Matrix(self.hight, other.width, lambda i, j: 0)

        for i in range(self.hight):
            for j in range(other.width):
                for k in range(self.width):
                    res.data[i][j] += self.data[i][k] * other.data[k][j]

        return res

    
    def __pow__(self, power):
        if power < -1:
            raise ValueError("Power of matrix should be >= 1")
        elif power == -1:
            return self.getInversed()

        res = copy.deepcopy(self)
        for i in range(power - 1):
            res = res * self

        return res


    def getTransposed(self):
        res = copy.deepcopy(self)

        for i in range(res.hight):
            for j in range(res.width):
                res.data[i][j], res.data[j][i] = res.data[j][i], res.data[i][j]
        return res

    
    def getInversed(self):
        matrix = copy.deepcopy(self)

        augmented_matrix = [
            [
                matrix.data[i][j] if j < matrix.hight else int(i == j - matrix.hight)
                for j in range(2 * matrix.hight)
            ]
            for i in range(matrix.hight)
        ]
        for i in range(matrix.hight):
            pivot = augmented_matrix[i][i]
            if pivot == 0:
                raise ValueError("Matrix is not invertible")
            for j in range(2 * matrix.hight):
                augmented_matrix[i][j] /= pivot
            for j in range(matrix.hight):
                if i != j:
                    scalar = augmented_matrix[j][i]
                    for k in range(2 * matrix.hight):
                        augmented_matrix[j][k] -= scalar * augmented_matrix[i][k]
        inverse = [
            [augmented_matrix[i][j] for j in range(matrix.hight, 2 * matrix.hight)]
            for i in range(matrix.hight)
        ]
        return Matrix(inverse)


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

# Return Lower triangular, Upper triangular and Diagonal matrices
def decomposeMatrix(matrix_in: Matrix):
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

    return L, U, D

    
def solveJacobi(matrix_in: Matrix, right_side: list, iter_amnt: int, answer: list = None):
    L, U, D = decomposeMatrix(matrix_in)
    R = (D ** (-1)) * (-1) * (L + U)
    F = (D ** (-1)) * right_side
    X = [1 for i in range(L.width)] 

    iterations = [i for i in range(1, iter_amnt + 1)]
    diff = []

    for i in range(iter_amnt):
        tmp = R * X
        for i in range(L.width):
            X[i] = tmp[i] + F[i]

        if answer != None:
            diff_tmp = 0
            for i in range(L.width):
                diff_tmp += (X[i] - answer[i]) ** 2
            diff.append(diff_tmp**0.5)

    if answer != None:
        plt.title("Jacobi method")
        plt.ylabel("Difference")
        plt.xlabel("Iterations")
        
        plt.minorticks_on()
        plt.grid(which='major')
        plt.grid(which='minor', linestyle=':')
        
        plt.plot(iterations, diff)
        plt.savefig("images/Jacobi.png")
        plt.show()
    
    return X


def solveSeidel(matrix_in: Matrix, right_side: list, iter_amnt: int, answer: list = None):
    L, U, D = decomposeMatrix(matrix_in)
    R = ((L + D) ** (-1)) * (-1) * U
    F = ((L + D) ** (-1)) * right_side
    X = [1 for i in range(L.width)]

    iterations = [i for i in range(1, iter_amnt + 1)]
    diff = []

    for i in range(iter_amnt):
        tmp = R * X
        for i in range(L.width):
            X[i] = tmp[i] + F[i]
        
        if answer != None:
            diff_tmp = 0
            for i in range(L.width):
                diff_tmp += (X[i] - answer[i]) ** 2
            diff.append(diff_tmp**0.5)

    if answer != None:
        plt.title("Siedel method")
        plt.ylabel("Difference")
        plt.xlabel("Iterations")
        
        plt.minorticks_on()
        plt.grid(which='major')
        plt.grid(which='minor', linestyle=':')

        plt.plot(iterations, diff)
        plt.savefig("images/Siedel.png")
        plt.show()

    return X

 
n = 10
matrix = Matrix(n, n, lambda i, j: 1 * (j == i) + (1 / (j + i)) * (i != j))
right_side = [1 / i for i in range(1, n + 1)]

resGauss = solveGauss(matrix, right_side)

resJacobi = solveJacobi(matrix, right_side, 20, resGauss)
resSeidel = solveSeidel(matrix, right_side, 20, resGauss)
