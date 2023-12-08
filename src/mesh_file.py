import numpy as np
from itertools import product

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

def total_nodes(num_points, nbounds, direction):
    """ This function calculates the total number of nodes in the grid.
    Parameters:
    num_points (list): A list of integers that specify the number of points in each dimension.
    nbounds (int): The number of boundary nodes.
    direction (int): The direction of the boundary.

    Returns:
    totnode (int): The total number of nodes in the grid.
    """

    totnode = 1
    points = num_points.copy()
    

    if direction == 1:
        points[0] = num_points[0] + nbounds
    
    elif direction == 2:
        points[1] = num_points[1] + nbounds

    elif direction == 3:
        points[2] = num_points[2] + nbounds

    for i in range(len(points)):
        totnode *= points[i]
    return totnode


def create_grid(total_nodes, num_points, dx, width, thick, length):
    """ This function specifies the locations of the material points.
    Parameters:
    total_nodes (int): The total number of nodes in the grid.
    num_points (list): A list of integers that specify the number of points in each dimension.
    dx (float): Discretization along x.
    width (float): The width of the material.
    thick (float): The thickness of the material.
    length (float): The length of the material.
    Returns:
    x (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    alflag (numpy.ndarray): A numpy array with dimensions (total_nodes, 1) that sets the total number of material pointsi
    in the family to 100.
    """

    alflag = np.zeros((total_nodes, 1), dtype=int)
    dimension = len(num_points)
    x = np.zeros((total_nodes, dimension), dtype=float)

    indices = product(range(num_points[2]), range(num_points[1]), range(num_points[0]))

    for nnum, (i, j, k) in enumerate(indices):
        val_x = (dx / 2.0) + k * dx
        val_y = -0.5 * width + (dx / 2.0) + j * dx
        val_z = -0.5 * thick + (dx / 2.0) + i * dx

        x[nnum, 0] = val_x
        x[nnum, 1] = val_y
        x[nnum, 2] = val_z

        if val_x > (length - dx):
            alflag[nnum, 0] = 1


    return x, alflag

def boundary_region(x, num_points, dx, width, thick):
    """ This function specifies the material points in the boundary region.
    Parameters:
    x (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    num_points (list): A list of integers that specify the number of points in each dimension.
    dx (float): Discretization along x.
    width (float): The width of the material.
    thick (float): The thickness of the material.
    length (float): The length of the material.
    Returns:
    x (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points
    and boundary.
    """

    ndivz, ndivy, nbnds = num_points[2], num_points[1], 3

    indices = product(range(1, ndivz + 1), range(1, ndivy + 1), range(1, nbnds + 1))

    for nnum, (i, j, k) in enumerate(indices, start=10000):
        coordx = - (dx / 2.0) - (k - 1) * dx
        coordy = -0.5 * width + (dx / 2.0) + (j - 1) * dx
        coordz = -0.5 * thick + (dx / 2.0) + (i - 1) * dx
        x[nnum, 0] = coordx
        x[nnum, 1] = coordy
        x[nnum, 2] = coordz
        nnum += 1

    return x

def horizon(total_nodes, u, delta):
    """ This function computes the horizon of each material point.
    Parameters:
    total_nodes (int): The total number of nodes in the grid.
    u (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    delta (float): The horizon of each material point x_k.
    Returns:
    numfam (numpy.ndarray): A numpy array with dimensions (total_nodes, 1) that specifies Number of family 
    members of each material point.
    pointfam (numpy.ndarray): Index array to find the family members in nodefam array
    nodefam (numpy.ndarray): A numpy array with dimensions (10000000, 1) that specifies the family members of each material point.
    """

    numfam = np.zeros((total_nodes, 1), dtype=int)
    pointfam = np.zeros((total_nodes, 1), dtype=int)
    nodefam = np.zeros((10000000, 1), dtype=int)

    for i in range(total_nodes):
        if i == 0:
            pointfam[i, 0] = 1
        else:
            pointfam[i, 0] = pointfam[i-1, 0] + numfam[i-1, 0]

        for j in range(total_nodes):
            idist = np.sqrt((u[i, 0] - u[j, 0])**2 + (u[i, 1] - u[j, 1])**2 + (u[i, 2] - u[j, 2])**2)
            if i != j and idist <= delta:
                numfam[i-1, 0] += 1
                nodefam[pointfam[i-1, 0] + numfam[i-1, 0] - 1, 0] = j

    return numfam, pointfam, nodefam

def surface_correction(total_nodes, u, numfam, modulus):

    sedload1 = 0.6 * modulus * 1.0e-6

    disp = np.zeros((total_nodes, 3), dtype=float)
    disp[:, 0] = 0.001 * u[:, 0]
    disp[:, 1:] = 0.0

    strain_dens = np.zeros((total_nodes, 3), dtype=float)
    fncst = np.zeros((total_nodes, 3), dtype=float)

    for i in range(total_nodes):
        for j in range(numfam[i, 0]):
            cnode = nodefam[pointfam[i, 0] + j-1, 0]
            idist = np.sqrt(np.sum((coord[cnode, :] - coord[i, :])**2))
            nlength = np.sqrt(np.sum((coord[cnode, :] + disp[cnode, :] - coord[i, :] - disp[i, :])**2))

            if idist <= delta - radij:
                fac = 1.0
            elif idist <= delta + radij:
                fac = (delta + radij - idist) / (2.0 * radij)
            else:
                fac = 0.0

            strain_dens[i, 0] += 0.5 * 0.5 * bc * ((nlength - idist) / idist)**2 * idist * vol * fac

        fncst[i, 0] = sedload1 / strain_dens[i, 0]
    
    return strain_dens, fncst
            
        






    
    


