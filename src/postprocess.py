import numpy as np
import matplotlib.pyplot as plt

def extract_values(coord, disp, y_value, z_value):
    """Extract the displacement values for a given y and z value.
    Parameters:
    -----------
    coord : array_like
        Array of coordinates (x, y, z).
    disp : array_like
        Array of displacements (ux, uy, uz).
    y_value : float
        Value of y to extract.
    z_value : float
        Value of z to extract.
    Returns:
    --------
    x : array_like
        Array of x coordinates.
    disp : array_like   
        Array of displacements (ux, uy, uz).
    """

    mask = (coord[:, 1] == y_value) & (coord[:, 2] == z_value)
    return coord[mask, 0], disp[mask]

def convert(s):
    return float(s.decode().replace('D', 'E'))

data = np.genfromtxt('horizontal_dispsbt.txt', delimiter=None, converters={i: convert for i in range(9)})
data_vert = np.genfromtxt('transverse_dispsbt.txt', delimiter=None, converters={i: convert for i in range(9)})
data_long = np.genfromtxt('vertical_dispsbt.txt', delimiter=None, converters={i: convert for i in range(9)})

data_beam = np.genfromtxt('coord_disp_pd_nt.txt', delimiter=None, converters={i: convert for i in range(3)})

print(data.shape)

# Split the data into separate arrays
grid = data_beam[:, 0]  # x, y, z
disp = data_beam[:, 1]  # ux, uy, uz
analytical_disp = data_beam[:, 2]  # ux, uy, uz
print(grid.shape)

plt.figure()
plt.xlabel('x (m)')
plt.ylabel('Displacement (m)')
plt.title('Transverse Displacement of the bar along y-axis')
plt.plot(grid, disp, label='Peridynamic Solution')
plt.plot(grid, analytical_disp, label='analytical displacemment')
plt.legend(loc='best')
plt.show()
