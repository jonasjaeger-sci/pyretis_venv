2D WCA RETIS simulation
=======================

Simulation
----------
task = retis
steps = 30
interfaces = [1.24, 1.34, 1.40, 1.46, 1.52, 1.54, 1.64, 1.74]

System
------
units = real

Engine settings
---------------
class = lammps
lmp = ~/lammps/build/lmp
input_path = lammps_input
subcycles = 1
lammps_format = data
extra_files = ['dw-wca.in']

TIS settings
------------
freq = 0.5
maxlength = 20000
aimless = True
allowmaxlength = False
zero_momentum = True
rescale_energy = 1
sigma_v = -1
seed = 0

Initial-path
------------
method = load
load_folder = load

RETIS settings
--------------
swapfreq = 0.5
relative_shoots = None
nullmoves = True
swapsimul = True

Output settings
---------------
pathensemble-file = 1
screen = 10
order-file = 1
energy-file = 1
trajectory-file = 1
