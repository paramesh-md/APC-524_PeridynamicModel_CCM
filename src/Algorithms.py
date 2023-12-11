import numpy as np
import math
import mesh_file as mf

def mass_vector(total_nodes, delta, bc, dx, **kwargs):

    mass_array = np.zeros((total_nodes, 3), dtype=int)

    mass_array[:, :3] = 0.25 * dt * dt * (4/3) * math.pi * (delta * delta * delta) * bc / dx

    return mass_array

def solver(total_nodes, alflag, dx, **kwargs):

    dtemp = kwargs['dtemp']
    pressure = kwargs['applied_pressure']

    bforce = np.zeros((total_nodes, 1), dtype=float)
    disp = np.zeros((total_nodes, 3), dtype=float)
    velocity = np.zeros((total_nodes, 3), dtype=float)
    

    bforce[alflag[:, 0] == 1, 0] = pressure / dx

    for t in rang(1, **kwargs['time_steps']+1):

        for n in range(total_nodes):
            pforce = np.zeros((total_nodes, 1), dtype=float)

            for j in range(1, numfam[n, 0]+1):
                cnode = nodefam[pointfam[n, 0] + j-1, 0]
                idist = mf.calculate_idist(u, n, cnode)
                nlength = mf.calculate_nlength(u, disp, n, cnode)

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
                        theta = math.atan(abs(u[cnode, 1] - u[i, 1]) / abs(u[cnode, 0] - u[i, 0]))

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
                    theta = math.atan(abs(u[cnode, 1] - u[i, 1]) / abs(u[cnode, 0] - u[i, 0]))
                    phi = math.acos(abs(u[cnode, 2] - u[i, 2]) / idist)

                    scx = (fncst[i, 0] + fncst[cnode, 0]) / 2.0
                    scy = (fncst[i, 1] + fncst[cnode, 1]) / 2.0
                    scz = (fncst[i, 2] + fncst[cnode, 2]) / 2.0
                    scr = 1.0 / (((math.cos(theta) * math.sin(phi))**2 / scx**2) + ((math.sin(theta) * math.sin(phi))**2 / scy**2) + (math.cos(phi)**2 / scz**2))
                    scr = math.sqrt(scr)

                # Calculation of the peridynamic force in x, y and z directions 
                # acting on a material point i due to a material point j
                dforce1 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (u[cnode, 0] + disp[cnode, 0] - u[i, 0] - disp[i, 0]) / nlength 
                dforce2 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (u[cnode, 1] + disp[cnode, 1] - u[i, 1] - disp[i, 1]) / nlength
                dforce3 = bc * ((nlength - idist) / idist - (alpha * dtemp)) * vol * scr * fac * (u[cnode, 2] + disp[cnode, 2] - u[i, 2] - disp[i, 2]) / nlength

                pforce[i, 0] = pforce[i, 0] + dforce1      
                pforce[i, 1] = pforce[i, 1] + dforce2  
                pforce[i, 2] = pforce[i, 2] + dforce3

                





    return bforce





