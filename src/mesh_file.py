import math
import numpy as np
from itertools import product

def mat_parameters(**kwargs):
    """ This function calculates the material parameters.
    Parameters:
    length (float): The length of the material.
    width (float): The width of the material.
    thick (float): The thickness of the material.
    num_points (list): A list of integers that specify the number of points in each dimension.
    modulus (float): The modulus of the material.
    Returns:
    dx (float): Discretization along x.
    dy (float): Discretization along y.
    dz (float): Discretization along z.
    delta (float): The horizon of each material point x_k.
    bc (float): The bond constant.
    area (float): The area of each material point.
    volume (float): The volume of each material point.
    """
    
    kwargs['dx'] = kwargs['length'] / kwargs['num_points'][0]
    kwargs['dy'] = kwargs['width'] / kwargs['num_points'][1]
    kwargs['dz'] = kwargs['thick'] / kwargs['num_points'][2]
    
    delta = 3.015 * kwargs['dx']
    bc = (12.0 * kwargs['modulus'])/(math.pi*(delta*delta*delta*delta))

    area = kwargs['dx'] * kwargs['dx']
    volume = area * kwargs['dx']

    kwargs['delta'] = delta
    kwargs['bond_constant'] = bc
    kwargs['area'] = area
    kwargs['volume'] = volume

    return kwargs

def total_nodes(**kwargs):
    """ This function calculates the total number of nodes in the grid.
    Parameters:
    num_points (list): A list of integers that specify the number of points in each dimension.
    nbounds (int): The number of boundary nodes.
    Returns:
    totnode (int): The total number of nodes in the grid.
    """

    totnode = 1
    points = kwargs['num_points'].copy()
    

    if kwargs['boundary_direction'] == 1:
        points[0] = kwargs['num_points'][0] + kwargs['nbounds']
    
    elif kwargs['boundary_direction'] == 2:
        points[1] = kwargs['num_points'][1] + kwargs['nbounds']

    elif kwargs['boundary_direction'] == 3:
        points[2] = kwargs['num_points'][2] + kwargs['nbounds']

    for i in range(len(points)):
        totnode *= points[i]
    return totnode


def create_grid(total_nodes, **kwargs):
    """ This function creates the grid of material points.
    Parameters:
    total_nodes (int): The total number of nodes in the grid.
    num_points (list): A list of integers that specify the number of points in each dimension.
    dx (float): Discretization along x.
    width (float): The width of the material.
    thick (float): The thickness of the material.
    length (float): The length of the material.
    Returns:
    x (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    alflag (numpy.ndarray): A numpy array with dimensions (total_nodes, 1) that specifies the boundary.
    """

    alflag = np.zeros((total_nodes, 1), dtype=int)
    dimension = len(kwargs['num_points'])
    x = np.zeros((total_nodes, dimension), dtype=float)

    dx = kwargs['dx']

    indices = product(range(kwargs['num_points'][2]), range(kwargs['num_points'][1]), range(kwargs['num_points'][0]))

    for nnum, (i, j, k) in enumerate(indices):

        val_x = (dx / 2.0) + k * dx
        val_y = -0.5 * kwargs['width'] + (dx / 2.0) + j * dx
        val_z = -0.5 * kwargs['thick'] + (dx / 2.0) + i * dx

        x[nnum, 0] = val_x
        x[nnum, 1] = val_y
        x[nnum, 2] = val_z

        if val_x > (kwargs['length'] - dx):
            alflag[nnum, 0] = 1


    return x, alflag

def boundary_region(x, **kwargs):
    """ This function creates the boundary region of material points.
    Parameters:
    x (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    num_points (list): A list of integers that specify the number of points in each dimension.
    dx (float): Discretization along x.
    width (float): The width of the material.
    thick (float): The thickness of the material.
    length (float): The length of the material.
    Returns:
    x (numpy.ndarray): An updated numpy array with dimensions (total_nodes, 3) that specifies the locations of the material and boundary points.
    """

    ndivz, ndivy, nbnds = kwargs['num_points'][2], kwargs['num_points'][1], 3

    indices = product(range(1, ndivz + 1), range(1, ndivy + 1), range(1, nbnds + 1))

    dx = kwargs['dx']

    for nnum, (i, j, k) in enumerate(indices, start=10000):
        coordx = - (dx / 2.0) - (k - 1) * dx
        coordy = -0.5 * kwargs['width'] + (dx / 2.0) + (j - 1) * dx
        coordz = -0.5 * kwargs['thick'] + (dx / 2.0) + (i - 1) * dx
        x[nnum, 0] = coordx
        x[nnum, 1] = coordy
        x[nnum, 2] = coordz
        nnum += 1

    return x

def calculate_idist(u, i, j):
    """ This function computes the Euclidean distance between two points u[i] and u[j].
    Parameters:
    u (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    i (int): The index of the first point.
    j (int): The index of the second point.
    Returns:
    idist (float): The Euclidean distance between u[i] and u[j].
    """

    idist = np.sqrt((u[i, 0] - u[j, 0])**2 + (u[i, 1] - u[j, 1])**2 + (u[i, 2] - u[j, 2])**2)
    return idist

def calculate_nlength(u, disp, i, j):
    """ This function computes the Euclidean distance between two points u[i] and u[j].
    Parameters:
    u (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    disp (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the displacements of the material points.
    i (int): The index of the first point.
    j (int): The index of the second point.
    Returns:
    nlength (float): The Euclidean distance between u[i] and u[j].
    """
    nlength = np.sqrt((u[i, 0] + disp[i, 0] - u[j, 0] - disp[j, 0])**2 + (u[i, 1] + disp[i, 1] - u[j, 1] - disp[j, 1])**2 + (u[i, 2] + disp[i, 2] - u[j, 2] - disp[j, 2])**2)
    return nlength


def horizon(total_nodes, u, **kwargs):
    """ This function determines the material points inside the horizon of each material point.
    Parameters:
    total_nodes (int): The total number of nodes in the grid.
    u (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    delta (float): The horizon of each material point x_k.
    Returns:
    numfam (numpy.ndarray): A numpy array with dimensions (total_nodes, 1) that specifies the number of material points inside the horizon of each material point.
    pointfam (numpy.ndarray): A numpy array with dimensions (total_nodes, 1) that specifies the starting index of the material points inside the horizon of each material point.
    nodefam (numpy.ndarray): A numpy array with dimensions (10000000, 1) that specifies the material points inside the horizon of each material point.
    """

    delta = kwargs['delta']

    numfam = np.zeros((total_nodes, 1), dtype=int)
    pointfam = np.zeros((total_nodes, 1), dtype=int)
    nodefam = np.zeros((10000000, 1), dtype=int)

    for i in range(total_nodes):
        if i == 0:
            pointfam[i, 0] = 1
        else:
            pointfam[i, 0] = pointfam[i-1, 0] + numfam[i-1, 0]

        for j in range(total_nodes):
            idist = calculate_idist(u, i, j)
            if i != j and idist <= delta:
                numfam[i, 0] += 1
                nodefam[pointfam[i, 0] + numfam[i, 0] - 1, 0] = j

    return [numfam, pointfam, nodefam]

def surface_correction(total_nodes, u, mat_family, dir, **kwargs):
    """ This function computes the surface correction factor for each material point.
    Parameters:
    total_nodes (int): The total number of nodes in the grid.
    u (numpy.ndarray): A numpy array with dimensions (total_nodes, 3) that specifies the locations of the material points.
    mat_family (list): A list of numpy arrays that specifies the material points inside the horizon of each material point.
    dir (int): The direction of loading.
    Returns:
    fncst (numpy.ndarray): A numpy array with dimensions (total_nodes, 1) that specifies the surface correction factor for each material point.
    """
    
    numfam, pointfam, nodefam = mat_family  # Unpack the arrays from horizon function

    dx = kwargs['dx']
    delta = kwargs['delta']
    bc = kwargs['bond_constant']
    volume = kwargs['volume']
    loading_direction = dir
    
    sedload1 = 0.6 * kwargs['modulus'] * 1.0e-6   # Strain energy density from CCM

    radij = 0.5 * delta
    disp = np.zeros((total_nodes, 3), dtype=float)
    
    if loading_direction == 0:
        disp[:, 0] = 0.001 * u[:, 0]
        disp[:, 1:] = 0.0

    elif loading_direction == 1:
        disp[:, 1] = 0.001 * u[:, 1]
        disp[:, 0] = 0.0
        disp[:, 2] = 0.0

    elif loading_direction == 2:
        disp[:, 2] = 0.001 * u[:, 2]
        disp[:, 0] = 0.0
        disp[:, 1] = 0.0
        
    strain_dens = np.zeros((total_nodes, 1), dtype=float)
    fncst = np.zeros((total_nodes, 1), dtype=float)

    for i in range(total_nodes):
        for j in range(1,numfam[i, 0]+1):
            cnode = nodefam[pointfam[i, 0] + j-1, 0]

            if i == cnode:
                print(i)
                continue

            idist = np.sqrt(np.sum((u[cnode, :] - u[i, :])**2))
            nlength = np.sqrt(np.sum((u[cnode, :] + disp[cnode, :] - u[i, :] - disp[i, :])**2))

            if idist <= delta - radij:
                fac = 1.0
            elif idist <= delta + radij:
                fac = (delta + radij - idist) / (2.0 * radij)
            else:
                fac = 0.0
            
            strain_dens[i, 0] += 0.5 * 0.5 * bc * ((nlength - idist) / idist)**2 * idist * volume * fac
        
        
        fncst[i, 0] = sedload1 / strain_dens[i, 0]
    
    return fncst
            
        






    
    


