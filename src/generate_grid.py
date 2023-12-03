import numpy as np
import matplotlib.pyplot as plt


def plot4d(data):
    """
    This function plots a 3D scatter plot of the input data. The color of the points is determined by the data values.
    The size of the points is determined by a mask that is True where the data is greater than 0.01 and False otherwise.
    The function can handle 2D or 3D data. If the data is 2D, a new array of zeros with the same shape as the data is created for the z-coordinates of the points.
    The plot is saved as a PNG image and then displayed.

    Parameters:
    data (numpy.ndarray): The input data to be plotted. It must be a 2D or 3D numpy array.

    Raises:
    ValueError: If the data is not 2D or 3D.
    """
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(projection="3d")
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    mask = data > 0.01
    idx = np.arange(int(np.prod(data.shape)))
    coordinates = np.unravel_index(idx, data.shape)
    if len(data.shape) == 2:
        x, y = coordinates
        z = np.zeros_like(x)
    elif len(data.shape) == 3:
        x, y, z = coordinates
    else:
        raise ValueError("Data must be 2D or 3D")
    ax.scatter(x, y, c=data.flatten(), s=10.0 * mask, edgecolor="face", alpha=0.2, marker="o", cmap="magma", linewidth=0)
    plt.tight_layout()
    plt.savefig("test_scatter_4d.png", dpi=250)
    plt.show()


def create_grid(bounds, num_points):
    """
    This function creates a grid of points in n dimensions. The grid is created using numpy's linspace and meshgrid functions.
    The bounds of the grid in each dimension are specified by the bounds parameter, which must contain 2, 4, or 6 elements.
    The number of points in each dimension is specified by the num_points parameter.

    Parameters:
    bounds (list): A list of 2, 4, or 6 elements that specify the bounds of the grid in each dimension.
    num_points (int): The number of points in each dimension.

    Returns:
    grid (list): A list of numpy arrays that represent the coordinates of the points in the grid.

    Raises:
    ValueError: If bounds does not contain 2, 4, or 6 elements.
    """

    dimensions = len(bounds) // 2
    if dimensions not in [1, 2, 3]:
        raise ValueError("bounds must contain 2, 4, or 6 elements")

    coordinates = [np.linspace(bounds[i*2], bounds[i*2+1], num_points) for i in range(dimensions)]
    grid = np.meshgrid(*coordinates, indexing='ij')

    return grid

        
        
#x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

num_points = 5
bounds = [0, 1, 0, 1]

grid = create_grid(bounds, num_points)
density_matrix = grid[0]**2 + grid[1]**2
print(density_matrix)
plot4d(density_matrix)
    

# if __name__ == "__main__":

#     assert np.all(x[:,0,0] == x_)
#     assert np.all(y[0,:,0] == y_)
#     assert np.all(z[0,0,:] == z_)


#     #X, Y, Z = np.meshgrid(x_, y_, z_, indexing="ij")
#     density_matrix = x**2 + y**2 + z**2
#     print(x_, y_, z_)
#     print(density_matrix)
#     plot4d(density_matrix)