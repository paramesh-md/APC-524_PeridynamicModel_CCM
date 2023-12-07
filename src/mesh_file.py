import numpy as np

def mat_parameters(length, width, thick, num_points):
    """ This function calculates the material parameters.
    Parameters:
    length (float): The length of the material.
    width (float): The width of the material.
    thick (float): The thickness of the material.
    num_points (list): A list of integers that specify the number of points in each dimension.
    Returns:
    dx (float): Discretization along x.
    dy (float): Discretization along y.
    dz (float): Discretization along z.
    delta (float): The horizon of each material point x_k.
    area (float): The area of x_k.
    volume (float): The volume of x_k.
    """
    
    dx = length / num_points[0]
    dy = width / num_points[1]
    dz = thick / num_points[2]
    
    delta = 3.015 * dx

    area = dx * dx
    volume = area * dx

    return [dx, dy, dz, delta, area, volume]

def total_nodes(num_points, nbounds):
    """ This function calculates the total number of nodes in the grid.
    Parameters:
    num_points (list): A list of integers that specify the number of points in each dimension.
    nbounds (int): The number of boundary nodes.
    Returns:
    totnode (int): The total number of nodes in the grid.
    """

    totnode = 1
    for i in range(len(num_points)):
        totnode *= num_points[i]
    return totnode + nbounds

def number_family(total_nodes):
    """ This function calculates the total number of nodes in the grid.
    Parameters:
    total_nodes (int): The total number of nodes in the grid.
    Returns:
    numfam (numpy.ndarray): An empty numpy array with dimensions (total_nodes, 1).
    pointfam (numpy.ndarray): An empty numpy array with dimensions (total_nodes, 1).   
    """
    numfam = np.zeros((total_nodes, 1), dtype=int)
    pointfam = np.zeros((total_nodes, 1), dtype=int)
    alflag = np.zeros((total_nodes, 1), dtype=int)

    return numfam, pointfam


def create_grid(total_nodes, num_points, dx, width, thick, length):

    alflag = np.zeros((total_nodes, 1), dtype=int)
    dimension = len(num_points)
    nnum = 0
    x = np.zeros((total_nodes, dimension), dtype=float)
    for i in range(num_points[2]):
        for j in range(num_points[1]):
            for k in range(num_points[0]):

                val_x = (dx / 2.0) + k * dx
                val_y = -0.5 * width + (dx / 2.0) + j * dx
                val_z = -0.5 * thick + (dx / 2.0) + i * dx

                x[nnum, 0] = val_x
                x[nnum, 1] = val_y
                x[nnum, 2] = val_z
                nnum += 1

                if val_x > (length-dx):
                    alflag[nnum, 0] = 1

    return x, alflag
    


