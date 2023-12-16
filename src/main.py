import numpy as np
from src import mesh_file as mf
from src import global_variables as gv
from src import Algorithms as alg
#import Algorithms as alg

gv.global_vars = mf.mat_parameters(**gv.global_vars)

print(gv.global_vars)
grid = mf.total_nodes(**gv.global_vars)

u, alflag = mf.create_grid(grid, **gv.global_vars)
#np.save('expected_u.npy', u)
#np.save('expected_alflag.npy', alflag)

#print(len(u), len(alflag))

u_updated = mf.boundary_region(u, **gv.global_vars)
#np.save('expected_u_updated.npy', u)
mat_family = mf.horizon(grid, u_updated, **gv.global_vars)

numfam, pointfam, nodefam = mat_family
#np.save('expected_numfam.npy', numfam)
##np.save('expected_pointfam.npy', pointfam)
#force = alg.solver(grid, alflag, dx, **gv.global_vars)
#fncst = np.zeros((grid, 3), dtype=float)

fncst[:, 0] = np.squeeze(mf.surface_correction(grid, u_updated, mat_family, dir = 0, **gv.global_vars))
fncst[:, 1] = np.squeeze(mf.surface_correction(grid, u_updated, mat_family, dir = 1, **gv.global_vars))
fncst[:, 2] = np.squeeze(mf.surface_correction(grid, u_updated, mat_family, dir = 2, **gv.global_vars))
#np.savetxt("strain_dens.txt", strain_dens, fmt='%.4f')
#np.savetxt("fncst.txt", fncst, fmt='%.4f')

displacement = alg.solver(grid, u_updated, alflag, mat_family, fncst, **gv.global_vars)
print(displacement)



#np.savetxt("myfile.txt", u, fmt='%.4f')
#np.save('expected_u.npy', u)
#np.save('expected_alflag.npy', alflag)
#np.savetxt('bforce.txt', force, fmt='%.4f')
