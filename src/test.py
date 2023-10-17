import netCDF4 as nc
import numpy as np
import h5py


#E0==sqrt(KTe0 ne/epsilon)

te= 50 #eV
ne= 1e19 #m^-3
epsilon= 8.854187817e-12 #F/m
e= 1.60217662e-19 #C
K= 1.38064852e-23 #J/K
E0= np.sqrt(K*te*ne*11600/epsilon)
debbye =np.sqrt(K*te*11600*epsilon/(e**2 * ne))
print(debbye)
