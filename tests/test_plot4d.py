import pytest
import numpy as np
from generate_grid import plot4d

def test_plot4d_2d():
    data = np.random.rand(10, 10)
    try:
        plot4d(data)
    except ValueError:
        pytest.fail("plot4d raised ValueError unexpectedly!")

def test_plot4d_3d():
    data = np.random.rand(10, 10, 10)
    try:
        plot4d(data)
    except ValueError:
        pytest.fail("plot4d raised ValueError unexpectedly!")

def test_plot4d_invalid_data():
    data = np.random.rand(10, 10, 10, 10)
    with pytest.raises(ValueError):
        plot4d(data)