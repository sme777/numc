from utils import *
from unittest import TestCase
import numpy as np

"""
For each operation, you should write tests to test  on matrices of different sizes.
Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""
class TestAdd(TestCase):
    def test_small_add(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 2, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_small_add1(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(10, 10, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(10, 10, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(50, 50, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(50, 50, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_add1(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(100, 100, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(100, 100, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_add2(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(500, 500, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(500, 500, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(1000, 1000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(1000, 1000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_add1(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(5000, 5000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(5000, 5000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_add2(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(10000, 10000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(10000, 10000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_add_iterative(self):
        speed_ups = []
        for i in range(15):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(5000, 5000, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(5000, 5000, seed=1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            speed_ups.append(speed_up)
            self.assertTrue(is_correct)
        speed_ups = np.array([speed_ups])
        print('The fastest was ' + str(np.max(speed_ups)))
        print('The slowest was ' + str(np.min(speed_ups)))
        print('The mean was ' + str(np.mean(speed_ups)))
    

class TestSub(TestCase):
    def test_small_sub(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 2, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_small_sub1(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(5, 5, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(5, 5, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_small_sub2(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(10, 10, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(10, 10, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_sub(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(50, 50, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(50, 50, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_sub1(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(100, 100, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(100, 100, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_sub1(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(500, 500, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(500, 500, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_sub(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(1000, 1000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(1000, 1000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_sub1(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(5000, 5000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(5000, 5000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_sub2(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(10000, 10000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(10000, 10000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_sub_iterative(self):
        speed_ups = []
        for i in range(15):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(5000, 5000, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(5000, 5000, seed=1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            speed_ups.append(speed_up)
            self.assertTrue(is_correct)
        speed_ups = np.array([speed_ups])
        print('The fastest was ' + str(np.max(speed_ups)))
        print('The slowest was ' + str(np.min(speed_ups)))
        print('The mean was ' + str(np.mean(speed_ups)))

class TestAbs(TestCase):
    def test_small_abs(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_small_abs1(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(5, 5, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_small_abs2(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 10, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_abs(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(50, 50, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_abs1(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(100, 100, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_abs2(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(500, 500, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_abs(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1000, 1000, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_abs1(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(5000, 5000, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_abs2(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10000, 10000, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_abs_iterative(self):
        speed_ups = []
        for i in range(15):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(5000, 5000, seed=0)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "abs")
            speed_ups.append(speed_up)
            self.assertTrue(is_correct)
        speed_ups = np.array([speed_ups])
        print('The fastest was ' + str(np.max(speed_ups)))
        print('The slowest was ' + str(np.min(speed_ups)))
        print('The mean was ' + str(np.mean(speed_ups)))
    

class TestNeg(TestCase):
    def test_small_neg(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_small_neg1(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(5, 5, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_small_neg2(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 10, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_neg(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(50, 50, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_neg1(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(100, 100, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_neg2(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(500, 500, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)


    def test_large_neg(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(1000, 1000, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_neg1(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(5000, 5000, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_neg2(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(10000, 10000, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_neg_iterative(self):
        speed_ups = []
        for i in range(15):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(5000, 5000, seed=0)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "neg")
            speed_ups.append(speed_up)
            self.assertTrue(is_correct)
        speed_ups = np.array([speed_ups])
        print('The fastest was ' + str(np.max(speed_ups)))
        print('The slowest was ' + str(np.min(speed_ups)))
        print('The mean was ' + str(np.mean(speed_ups)))
    

class TestMul(TestCase):
    def test_small_mul(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 3, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(3, 4, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
        
    def test_medium_mul(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(105, 205, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(205, 300, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_mul(self):
        speed_ups = []
        for i in range(5):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(1005, 1005, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(1005, 1111, seed=1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
            self.assertTrue(is_correct)
            speed_ups.append(speed_up)
        speed_ups = np.array([speed_ups])
        print('The fastest was ' + str(np.max(speed_ups)))
        print('The slowest was ' + str(np.min(speed_ups)))
        print('The mean was ' + str(np.mean(speed_ups)))

    # def test_large_mul1(self):
    #     dp_mat1, nc_mat1 = rand_dp_nc_matrix(3000, 5000, seed=0)
    #     dp_mat2, nc_mat2 = rand_dp_nc_matrix(5000, 2000, seed=1)
    #     is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    # def test_large_mul2(self):
    #     dp_mat1, nc_mat1 = rand_dp_nc_matrix(10000, 8000, seed=0)
    #     dp_mat2, nc_mat2 = rand_dp_nc_matrix(8000, 11000, seed=1)
    #     is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)
    

class TestPow(TestCase):
    def test_small_pow(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        is_correct, speed_up = compute([dp_mat, 3], [nc_mat, 3], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_pow(self):
        # TODO: YOUR CODE HERE
        speed_ups = []
        for i in range(5):
            dp_mat, nc_mat = rand_dp_nc_matrix(100, 100, seed=0)
            print("I work")
            is_correct, speed_up = compute([dp_mat, 1000], [nc_mat, 1000], "pow")
            self.assertTrue(is_correct)
            speed_ups.append(speed_up)
        speed_ups = np.array(speed_ups)
        print('The fastest was ' + str(np.max(speed_ups)))
        print('The slowest was ' + str(np.min(speed_ups)))
        print('The mean was ' + str(np.mean(speed_ups)))

        
        
    def test_large_pow(self):
        # TODO: YOUR CODE HERE
        pass

class TestGet(TestCase):
    def test_get(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        rand_row = np.random.randint(dp_mat.shape[0])
        rand_col = np.random.randint(dp_mat.shape[1])
        self.assertEqual(round(dp_mat[rand_row][rand_col], decimal_places),
            round(nc_mat[rand_row][rand_col], decimal_places))

class TestSet(TestCase):
    def test_set(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        rand_row = np.random.randint(dp_mat.shape[0])
        rand_col = np.random.randint(dp_mat.shape[1])
        self.assertEquals(round(dp_mat[rand_row][rand_col], decimal_places),
            round(nc_mat[rand_row][rand_col], decimal_places))

class TestShape(TestCase):
    def test_shape(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        self.assertTrue(dp_mat.shape == nc_mat.shape)


class TestSlicing(TestCase):
    def test_set(self):
        b = nc.Matrix([[1, 2, 3], [4, 5, 6], [7,8,9]])
       # a[0]
        b[1]
