import sys
import time
import numpy as np 
sys.path.append('/Users/wangdong/Documents/research/diffusion_kinetics_code_2.0')
from one_dimensional_grid import UniformGrid
from one_dimensional_diffusion import ExplicitDiffusion
from chemistry import Chemistry 
from snapshot_processing import SnapshotProcessing 

"""
input
"""
# mesh:
PT_file_location = './'
PT_file_name = 'PT_Jupiter.csv'
N = 60
PG_min = 30.e+5
PG_max = 2000.e+5
# diffusion:
K = 1.e4*np.ones(N)
# chemistry:
mechanism_path = '/Users/wangdong/Documents/Cornell/2014_summer/cantera_input_files/my_mechanism.cti'
# jupiter, 10 times solar water
composition = 'H2:1.0,He:0.157,CH4:2.37e-3,H2O:9.8e-3,NH3:6.64e-4'
atol = 1e-15
rtol = 1e-9
# simulation
dt = 2000
T_step = 10
data_frequency = 100
restart = 1
# post-processing
fig_name = 'mole_fraction.png' 
species_names = ['H2O','CH4','C2H6','CO','CO2']
output_file_name = 'mole_fractions.dat'
header_variable_names = 'Depth(m), Pressure(Pa), Temperature(K), X_H2O, X_CH4, X_C2H6, X_CO, X_CO2'

"""
One dimensional mesh
"""
print "start generating mesh..."
# read in the referene data
PG_ref, T_ref, z_ref, rho_ref = np.loadtxt(PT_file_location+PT_file_name, comments='#', delimiter=",", skiprows=0, unpack=True)
# generate unform mesh
uniform_grid = UniformGrid(N,PG_min,PG_max,PG_ref, T_ref, -z_ref, rho_ref)
print "mesh generated.","N =",N


"""
Diffusion object
"""
print "Diffusion solver setup..."
diffusion_solver = ExplicitDiffusion(K,uniform_grid)

"""
chemsitry object
"""
print "chemistry solver setup..."
chemistry_obj = Chemistry(mechanism_path,composition,uniform_grid,atol,rtol)
print "Number of species:", chemistry_obj.gas.n_species, "Number of reactions:",chemistry_obj.gas.n_reactions
print "The absolute error is", chemistry_obj.network.atol, "The relative error is", chemistry_obj.network.rtol

"""
time evolution
"""
print "Initialize the composition..."
if restart == 0:
	# a new run, start from equilibrium abundances
    y = chemistry_obj.equilibrate()
    # write equilibrium values into files
    np.savetxt('./y_equilibrium.txt',y)
    t_0 = 0
else:
    # restart
    # read the starting time from file
    f=open('./time.txt','r')
    t_0=int(f.read())
    f.close()
    # read the array from file
    y = np.loadtxt('./y_t.txt')

print "start time evolution..."
t_start = time.time()
n_species = chemistry_obj.gas.n_species
for j in range(T_step):
    # diffusion for dt/2
    y = diffusion_solver.advance(y,n_species,dt/2)
    # chemical evolution for dt
    chemistry_obj.advance(y,dt)
    # diffusion for another dt/2
    y = diffusion_solver.advance(y,n_species,dt/2)
    print j
    # record data for every 100 time steps and the end of the simulation
    if (j+1)%data_frequency == 0 or j+1 == T_step:
        # save data for plot or restart
        t = t_0 + (j+1)*dt
        # save the time to file
        f = open('./time.txt', 'w')
        f.write(str(t))
        f.close()
        np.savetxt('./y_t.txt',y)
        print "snapshot saved at t =",t

print "end of time evolution"
print "The time it takes is ", time.time()-t_start, "s"

"""
post processing
"""
post_obj = SnapshotProcessing(y,chemistry_obj,uniform_grid)
post_obj.Plot(species_names,fig_name)
post_obj.SaveData(species_names,output_file_name,header_variable_names)






