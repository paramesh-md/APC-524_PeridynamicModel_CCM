import pytest
import numpy as np
from generate_grid import create_grid

def test_create_grid_1d():
    bounds = [0, 1]
    num_points = 5
    grid = create_grid(bounds, num_points)
    assert len(grid) == 1
    assert np.all(grid[0] == np.linspace(0, 1, num_points))

def test_create_grid_2d():
    bounds = [0, 1, 0, 1]
    num_points = 5
    grid = create_grid(bounds, num_points)
    assert len(grid) == 2
    assert np.all(grid[0][:,0] == np.linspace(0, 1, num_points))
    assert np.all(grid[1][0,:] == np.linspace(0, 1, num_points))

def test_create_grid_3d():
    bounds = [0, 1, 0, 1, 0, 1]
    num_points = 5
    grid = create_grid(bounds, num_points)
    assert len(grid) == 3
    assert np.all(grid[0][:,0,0] == np.linspace(0, 1, num_points))
    assert np.all(grid[1][0,:,0] == np.linspace(0, 1, num_points))
    assert np.all(grid[2][0,0,:] == np.linspace(0, 1, num_points))

def test_create_grid_invalid_bounds():
    bounds = [0, 1, 0]
    num_points = 5
    with pytest.raises(ValueError):
        create_grid(bounds, num_points)