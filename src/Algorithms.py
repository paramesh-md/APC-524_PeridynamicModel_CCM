import numpy as np
import warnings
import math
from src import mesh_file as mf

def mass_vector(total_nodes, **kwargs):
    """This function calculates the mass vector for each material point in the domain. 
    The mass vector is calculated using the following equation:
    mass_vector = 0.25 * dt * dt * (4/3) * math.pi * (delta * delta * delta) * bond_constant / dx
    Parameters:
    total_nodes (int) : Total number of nodes in the domain
    **kwargs : Keyword arguments passed from the main function
    Returns:
    mass_array (numpy array) : Mass vector for each material point in the domain
    """

    delta = kwargs['delta']
    bc = kwargs['bond_constant']
    dx = kwargs['dx']

    dt = kwargs['total_time']/kwargs['time_steps']

    mass_array = np.zeros((total_nodes, 3), dtype=int)

    mass_array[:, :3] = 0.25 * dt * dt * (4/3) * math.pi * (delta * delta * delta) * bc / dx

    return mass_array

def solver(total_nodes, u, alflag, mat_family, fncst, **kwargs):
    """This function solves the equations of motion using the adaptive dynamic relaxation method.
    Parameters:
    total_nodes (int) : Total number of nodes in the domain
    u (numpy array) : Material points in the model
    alflag (numpy array) : Flag for each material point
    mat_family (tuple) : Tuple containing the number of material points in the horizon of each material point, 
    the starting index of the material points in the horizon of each material point and the material points in the horizon of each material point
    fncst (numpy array) : Surface correction factor for each material point
    **kwargs : Keyword arguments passed from the main function
    Returns:
    disp (numpy array) : Displacement of each material point in the domain
    """


    dtemp = kwargs['dtemp']
    pressure = kwargs['applied_pressure']
    dx = kwargs['dx']
    vol = kwargs['volume']
    bc = kwargs['bond_constant']
    dt = kwargs['time_steps']
    delta = kwargs['delta']

    massvec = mass_vector(total_nodes, **kwargs)

    radij = 0.5 * delta
    alpha = 23.0e-6

    numfam, pointfam, nodefam = mat_family

    bforce = np.zeros((total_nodes, 3), dtype=float)
    disp = np.zeros((total_nodes, 3), dtype=float)

    velocity = np.zeros((total_nodes, 3), dtype=float)

    vel_half = np.zeros((total_nodes, 3), dtype=float)
    vel_half_old = np.zeros((total_nodes,3), dtype=float)

    pforce = np.zeros((total_nodes, 3), dtype=float)
    pforceold = np.zeros((total_nodes, 3), dtype=float)

    bforce[alflag[:, 0] == 1, 0] = pressure/dx

    for t in range(1, kwargs['time_steps']+1):

        for i in range(total_nodes):
            pforce = np.zeros((total_nodes, 3), dtype=float)

            for j in range(1, numfam[i, 0]+1):
                cnode = nodefam[pointfam[i, 0] + j-1, 0]
                idist = mf.calculate_idist(u, i, cnode)
                nlength = mf.calculate_nlength(u, disp, i, cnode)

                if idist <= delta - radij:
                    fac = 1.0
                
                elif idist <= delta + radij:
                    fac = (delta + radij - idist) / (2.0 * radij)
                else:
                    fac = 0.0
                
                if abs(u[cnode, 2] - u[i, 2]) <= 1e-10:

                    if abs(u[cnode, 1] - u[i, 1]) <= 1e-10:
                        theta = 0.0


                    elif abs(u[cnode, 0] - u[i, 0]) <= 1e-10:
                        theta = 90.0 * math.pi / 180.0

                    else:
                        theta = math.atan(abs(u[cnode, 1] - u[i, 1])/abs(u[cnode, 0] - u[i, 0]))
                        
                    phi = 90.0 * math.pi / 180.0

                    scx = (fncst[i,0] + fncst[cnode,0])/2
                    scy = (fncst[i,1] + fncst[cnode,1])/2
                    scz = (fncst[i,2] + fncst[cnode,2])/2

                    scr = 1.0 / (((math.cos(theta) * math.sin(phi))**2 / scx**2) + ((math.sin(theta) * math.sin(phi))**2 / scy**2) + (math.cos(phi)**2 / scz**2))
                    scr = math.sqrt(scr)
                
                elif abs(u[cnode, 0] - u[i, 0]) <= 1e-10 and abs(u[cnode, 1] - u[i, 1]) <= 1e-10:
                    scz = (fncst[i, 2] + fncst[cnode, 2]) / 2.0
                    scr = scz

                else:
                    if abs(u[cnode, 0] - u[i, 0]) <= 1e-10:
                        theta = 90.0 * math.pi / 180.0
                    else:
                        theta = math.atan(abs(u[cnode, 1] - u[i, 1]) / abs(u[cnode, 0] - u[i, 0]))
                    phi = math.acos(abs(u[cnode, 2] - u[i, 2]) / idist)

                    scx = (fncst[i, 0] + fncst[cnode, 0]) / 2.0
                    scy = (fncst[i, 1] + fncst[cnode, 1]) / 2.0
                    scz = (fncst[i, 2] + fncst[cnode, 2]) / 2.0
    
                    scr = 1.0 / (((math.cos(theta) * math.sin(phi))**2 / scx**2) + ((math.sin(theta) * math.sin(phi))**2 / scy**2) + (math.cos(phi)**2 / scz**2))
                    scr = math.sqrt(scr)

                # Calculation of the peridynamic force in x, y and z directions acting on a material point i due to a material point j
                dforce1 = bc*((nlength - idist)/idist - (alpha * dtemp))*vol*scr*fac*(u[cnode, 0] + disp[cnode, 0] - u[i, 0] - disp[i, 0])/nlength 
                dforce2 = bc*((nlength - idist)/idist - (alpha * dtemp))*vol*scr*fac*(u[cnode, 1] + disp[cnode, 1] - u[i, 1] - disp[i, 1])/nlength
                dforce3 = bc*((nlength - idist)/idist - (alpha * dtemp))*vol*scr*fac*(u[cnode, 2] + disp[cnode, 2] - u[i, 2] - disp[i, 2])/nlength

                pforce[i, 0] = pforce[i, 0] + dforce1      
                pforce[i, 1] = pforce[i, 1] + dforce2  
                pforce[i, 2] = pforce[i, 2] + dforce3

        # Adaptive Dynamic Relaxation

        cn = 0.0
        cn1 = 0.0
        cn2 = 0.0

        for i in range(total_nodes):

            if vel_half_old[i, 0] != 0.0:
                cn1 -= disp[i, 0] * disp[i, 0] * (pforce[i, 0]/massvec[i, 0] - pforceold[i, 0]/massvec[i, 0])/(dt*vel_half_old[i, 0])

            if vel_half_old[i, 1] != 0.0:
                cn1 -= disp[i, 1] * disp[i, 1] * (pforce[i, 1]/massvec[i, 1] - pforceold[i, 1]/massvec[i, 1])/(dt*vel_half_old[i, 1])

            if vel_half_old[i, 2] != 0.0:
                cn1 -= disp[i, 2] * disp[i, 2] * (pforce[i, 2]/massvec[i, 2] - pforceold[i, 2]/massvec[i, 2])/(dt*vel_half_old[i, 2])
                
            cn2 += disp[i, 0] * disp[i, 0]
            cn2 += disp[i, 1] * disp[i, 1]
            cn2 += disp[i, 2] * disp[i, 2]

        if cn2 != 0.0:
            if (cn1/cn2) > 0.0:
                cn = 2.0*math.sqrt(cn1/cn2)
            else:
                cn = 0.0
        else:
            cn = 0.0

        if cn > 2.0:
            cn = 1.9

        for val in range(total_nodes):

            if t == 1:
                vel_half[i, 0] = (dt/massvec[i, 0])*(pforce[i, 0] + bforce[i, 0])/2.0
                vel_half[i, 1] = (dt/massvec[i, 1])*(pforce[i, 1] + bforce[i, 1])/2.0
                vel_half[i, 2] = (dt/massvec[i, 2])*(pforce[i, 2] + bforce[i, 2])/2.0

            else:
                vel_half[i, 0] = ((2.0 - cn*dt)*vel_half_old[i, 0] + (2*dt/massvec[i, 0])*(pforce[i, 0] + bforce[i, 0]))/(2.0 + cn*dt)
                vel_half[i, 1] = ((2.0 - cn*dt)*vel_half_old[i, 1] + (2*dt/massvec[i, 1])*(pforce[i, 1] + bforce[i, 1]))/(2.0 + cn*dt)
                vel_half[i, 2] = ((2.0 - cn*dt)*vel_half_old[i, 2] + (2*dt/massvec[i, 2])*(pforce[i, 2] + bforce[i, 2]))/(2.0 + cn*dt)

            velocity[i, 0] = 0.5*(vel_half_old[i, 0] + vel_half[i, 0])
            velocity[i, 1] = 0.5*(vel_half_old[i, 1] + vel_half[i, 1])
            velocity[i, 2] = 0.5*(vel_half_old[i, 2] + vel_half[i, 2])

            disp[i, 0] = disp[i, 0] + vel_half[i, 0] * dt
            disp[i, 1] = disp[i, 1] + vel_half[i, 1] * dt
            disp[i, 2] = disp[i, 2] + vel_half[i, 2] * dt

            vel_half_old[i, 0] = vel_half[i, 0]
            vel_half_old[i, 1] = vel_half[i, 1]
            vel_half_old[i, 2] = vel_half[i, 2]

            pforceold[i, 0] = pforce[i, 0]
            pforceold[i, 1] = pforce[i, 1]
            pforceold[i, 2] = pforce[i, 2]
        

    return disp





