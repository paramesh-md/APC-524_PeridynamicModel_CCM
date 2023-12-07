import unittest
import numpy as np
import mesh_file as mf

class TestMeshFile(unittest.TestCase):
    def test_mat_parameters(self):
        length = 1.0
        width = 0.1
        thick = 0.1
        num_points = [100, 10, 10]
        dx, dy, dz, delta, area, volume = mf.mat_parameters(length, width, thick, num_points)
        self.assertEqual(dx, 0.01)
        # Add more assertions here

    def test_total_nodes(self):
        num_points = [103, 10, 10]
        nbounds = 3
        grid = mf.total_nodes(num_points, nbounds)
        self.assertEqual(grid, 10303)
        # Add more assertions here

    def test_number_family(self):
        grid = 10300
        numfam, pointfam = mf.number_family(grid)
        # Add assertions here

    def test_create_grid(self):
        num_points = [100, 10, 10]
        length = 1.0
        width = 0.1
        thick = 0.1
        dx = 0.01
        grid = 10300
        u, alflag = mf.create_grid(grid, num_points, dx, width, thick, length)
        self.assertEqual(u.shape, (grid, 3))
        # Add more assertions here

    def test_regression(self):
        num_points = [100, 10, 10]
        time_steps = 4000
        bounds = [0, 1, 0, 0.1, 0, 0.1]
        nbounds= 3

        length = 1.0
        width = 0.1
        thick = 0.1

        dx, dy, dz, delta, area, volume = mf.mat_parameters(length, width, thick, num_points)
        grid = mf.total_nodes(num_points, nbounds)
        numfam, pointfam = mf.number_family(grid)

        u, alflag = mf.create_grid(grid, num_points, dx, width, thick, length)

        # Compare u with expected result
        expected_u = np.load('expected_u.npy')
        np.testing.assert_array_almost_equal(u, expected_u, decimal=4)

        # Compare alflag with expected result
        expected_alflag = np.load('expected_alflag.npy')
        np.testing.assert_array_almost_equal(alflag, expected_alflag, decimal=4)

if __name__ == '__main__':
    unittest.main()