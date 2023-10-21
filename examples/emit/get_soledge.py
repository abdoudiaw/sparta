import numpy as np
import matplotlib.pyplot as plt

def get_sol_data():
        sys.path.insert(0, '/Users/42d/syndiag/WEST/NN_SOLEDGE/')
        import nn_learner
        from SOLEDEG3x import Soledge
        
        data_path = '/Users/42d/syndiag/WEST/NN_SOLEDGE/data_soledge.pt'
        soledge = Soledge(data_path)
        
        data = h5py.File("west_54695_totally_works.h5")
        z_m = data['z'][:]
        r_m = data['r'][:]

        te = [soledge.get_soledge_data(r_values, z_values).TEMP_E for r_values, z_values in zip(r_m, z_m)]
        ne = [soledge.get_soledge_data(r_values, z_values).DENS_E for r_values, z_values in zip(r_m, z_m)]
        
        sys.path.pop(0)
        if 'nn_learner' in sys.modules:
            del sys.modules['nn_learner']

        return np.asarray(te), np.asarray(ne)


