"""
'Global variables for the program'
'num_points: number of points in each direction'
'total_time: total time of simulation'
'time_steps: number of time steps'
'bounds: Model dimensions and coordinates of the boundary'
'nbounds: Number of boundary conditions'
'length: Length of the model'
'width: Width of the model'
'thick: Thickness of the model'
'modulus: Elastic modulus of the material'
'boundary_direction: Identify boundary points along x, y or z direction'
'loading_direction: Direction of loading'
'applied_pressure: Applied pressure'
'dtemp: Change in temperature'
'volume: Volume of the material point - Is calculated in the program'
'area: Area of the material point - Is calculated in the program'
'delta: Horizon of the material point - Is calculated in the program'
'dx: Grid spacing - Is calculated in the program'
'bond_constant: Bond constant - Is calculated in the program'
"""
global_vars = {
    'num_points':[100, 10, 10],
    'total_time':4000,
    'time_steps':4000,
    'bounds':[0, 1, 0, 0.1, 0, 0.1],
    'nbounds':3,
    'length':1.0,
    'width':0.1,
    'thick':0.1,
    'modulus':200.0,
    'boundary_direction':1,
    'loading_direction':[1, 2, 3],
    'applied_pressure':200.0,
    'dtemp':0.0,
    'volume':None,
    'area':None,
    'delta':None,
    'dx':None,
    'bond_constant':None}