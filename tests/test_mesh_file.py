import unittest
import numpy as np
import mesh_file

class TestMeshFile(unittest.TestCase):

    def test_mat_parameters(self):
        kwargs = {
            'length': 1,
            'width': 0.1,
            'thick': 0.1,
            'num_points': [100, 10, 10],
            'modulus': 200
        }
        result = mesh_file.mat_parameters(**kwargs)
        expected_dx = kwargs['length'] / kwargs['num_points'][0]
        expected_dy = kwargs['width'] / kwargs['num_points'][1]
        expected_dz = kwargs['thick'] / kwargs['num_points'][2]
        expected_delta = 3.015 * expected_dx
        expected_bc = (12.0 * kwargs['modulus'])/(np.pi*(expected_delta**4))
        self.assertEqual(result['dx'], 0.01)
        self.assertEqual(result['dy'], 0.01)
        self.assertEqual(result['dz'], 0.01)
        self.assertEqual(round(result['delta'],5), 0.03015)
        self.assertEqual(round(result['bond_constant'],4), 924511043.1546)

    def test_mat_parameters_regression(self):  

        kwargs = {
            'length': 1,
            'width': 0.1,
            'thick': 0.1,
            'num_points': [100, 10, 10],
            'modulus': 200.0
        }
        result = mesh_file.mat_parameters(**kwargs)

        expected_result = {
            'dx': 0.01,
            'dy': 0.01,
            'dz': 0.01,
            'delta': 0.03015,
            'bond_constant': 924511043.1546
        }
        for key in expected_result.keys():
            self.assertAlmostEqual(result[key], expected_result[key], places=4)

    def test_total_nodes(self):
        kwargs = {
            'num_points': [100, 10, 10],
            'nbounds':3,
            'boundary_direction': 1
        }
        result = mesh_file.total_nodes(**kwargs)
        expected_result = (kwargs['num_points'][0] + kwargs['nbounds'])*kwargs['num_points'][1]*kwargs['num_points'][2]
        self.assertEqual(result, expected_result)

    def test_total_nodes_regression(self):
        kwargs = {
            'num_points': [100, 10, 10],
            'nbounds':3,
            'boundary_direction': 1
        }
        result = mesh_file.total_nodes(**kwargs)
        # This is the expected result based on known inputs
        expected_result = 10300
        self.assertEqual(result, expected_result)

    def test_create_grid(self):
        num_points = [10, 5, 2]
        dx = 1.0
        dy = 2.0
        dz = 3.0
        nbounds = 1
        result = mesh_file.create_grid(num_points, dx, dy, dz, nbounds)
        expected_result = np.mgrid[-dx:dx*(num_points[0]+1), -dy:dy*(num_points[1]+1), -dz:dz*(num_points[2]+1)]
        np.testing.assert_array_equal(result, expected_result)

    def test_create_grid(self):
        total_nodes = 1000
        kwargs = {
            'num_points': [10, 10, 10],
            'dx': 1.0,
            'width': 10.0,
            'thick': 10.0,
            'length': 10.0
        }
        x, alflag = mesh_file.create_grid(total_nodes, **kwargs)
        # Check that x and alflag have the correct shape
        self.assertEqual(x.shape, (total_nodes, len(kwargs['num_points'])))
        self.assertEqual(alflag.shape, (total_nodes, 1))
        # Check that the values in x are correct
        for i in range(total_nodes):
            self.assertEqual(x[i, 0], (kwargs['dx'] / 2.0) + (i % kwargs['num_points'][0]) * kwargs['dx'])
            self.assertEqual(x[i, 1], -0.5 * kwargs['width'] + (kwargs['dx'] / 2.0) + ((i // kwargs['num_points'][0]) % kwargs['num_points'][1]) * kwargs['dx'])
            self.assertEqual(x[i, 2], -0.5 * kwargs['thick'] + (kwargs['dx'] / 2.0) + (i // (kwargs['num_points'][0] * kwargs['num_points'][1])) * kwargs['dx'])
        # Check that the values in alflag are correct
        for i in range(total_nodes):
            if x[i, 0] > (kwargs['length'] - kwargs['dx']):
                self.assertEqual(alflag[i, 0], 1)
            else:
                self.assertEqual(alflag[i, 0], 0)

    def test_create_grid_regression(self):
        total_nodes = 10300
        kwargs = {
            'num_points': [100, 10, 10],
            'nbounds': 3,
            'dx': 0.01,
            'width': 0.1,
            'thick': 0.1,
            'length': 1
        }
        x, alflag = mesh_file.create_grid(total_nodes, **kwargs)
        # Load expected results from a file or define them directly in the test
        expected_x = np.load('tests/expected_u.npy')
        #expected_alflag = np.load('expected_alflag.npy')
        np.testing.assert_array_equal(x, expected_x)
        #np.testing.assert_array_equal(alflag, expected_alflag)

    def test_boundary_region(self):
        total_nodes = 10300
        
        kwargs = {
            'num_points': [100, 10, 10],
            'nbounds': 3,
            'dx': 0.01,
            'width': 1,
            'thick': 0.1,
            'length': 0.1
        }
        x, alflag = mesh_file.create_grid(total_nodes, **kwargs)
        result = mesh_file.boundary_region(x, **kwargs)
        # Check that the values in result are correct
        for i in range(10000, result.shape[0]):
            self.assertEqual(result[i, 0], - (kwargs['dx'] / 2.0) - ((i - 10000) % 3 + 1 - 1) * kwargs['dx'])
            self.assertEqual(result[i, 1], -0.5 * kwargs['width'] + (kwargs['dx'] / 2.0) + (((i - 10000) // 3) % kwargs['num_points'][1] + 1 - 1) * kwargs['dx'])
            self.assertEqual(result[i, 2], -0.5 * kwargs['thick'] + (kwargs['dx'] / 2.0) + ((i - 10000) // (3 * kwargs['num_points'][1]) + 1 - 1) * kwargs['dx'])
    
    def test_calculate_idist(self):

        u = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])
        i = 1
        j = 2
        result = mesh_file.calculate_idist(u, i, j)
        expected_result = np.sqrt((u[i, 0] - u[j, 0])**2 + (u[i, 1] - u[j, 1])**2 + (u[i, 2] - u[j, 2])**2)
        self.assertEqual(result, expected_result)

    def test_calculate_nlength(self):
        u = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]])
        disp = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        i = 1
        j = 2
        result = mesh_file.calculate_nlength(u, disp, i, j)
        expected_result = np.sqrt((u[i, 0] + disp[i, 0] - u[j, 0] - disp[j, 0])**2 + (u[i, 1] + disp[i, 1] - u[j, 1] - disp[j, 1])**2 + (u[i, 2] + disp[i, 2] - u[j, 2] - disp[j, 2])**2)
        self.assertEqual(result, expected_result)
  

if __name__ == '__main__':
    unittest.main()
