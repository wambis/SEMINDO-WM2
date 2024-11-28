
#
# Imports 
import numpy as np
from io import StringIO as StringIO
#
# ---
#

class Parameter:
    def __init__(self, config):
        self.num_sspl, self.num_bspl, self.desc = config
        self.values = []
        self.grid = []
    def set_values(self, values):
        self.values = values
    def get_values(self):
        return self.values
    def set_grid(self, grid):
        self.grid = grid
    def get_grid(self):
        return self.grid
    def get_desc(self):
        return self.desc
        
class Discontinuity:
    def __init__(self, config):
        self.num_sspl, self.desc = config
        self.values = []
        self.grid = []
    def set_values(self, values):
        self.values = values
    def set_grid(self, grid):
        self.grid = grid

class ModelA3d:
    """
    A python class for manipulating UCB A3d (spherical spline + vertical b-spline) models
    """

    def __init__(self, model_file_name = []):
        # set the model file name
        self.file_name = model_file_name

        # set the populated flag to False
        self.populated = False

    def set_file_name(self, model_file_name):
        # set the model file name
        self.file_name = model_file_name

    def get_file_name(self):
        # return the model file name
        return self.file_name

    def set_bspl_knots(self, bspl_knots):
        self.bspl_knots = bspl_knots

    def get_bspl_knots(self):
        return self.bspl_knots

    def get_parameter_by_name(self, target):
        for p in self.parameters:
            if p.get_desc() == target:
                return p
        return None

    def set_parameters(self, parameter_structs):
        self.parameters = []
        # loop over parameter configuration / value structures
        for p_struct in parameter_structs:
            p_config, p_values = p_struct
            # initialize new Parameter object
            p_new = Parameter(p_config)
            p_new.set_values(p_values)
            self.parameters.append(p_new)
        self.num_parameters = len(parameter_structs)
    
    def set_discons(self, discon_structs):
        self.discons = []
        # loop over discontinuity configuration / value structures
        for d_struct in discon_structs:
            d_config, d_values = d_struct
            # initialize new Discontinuity object
            d_new = Discontinuity(d_config)
            d_new.set_values(d_values)
            self.discons.append(d_new)
        self.num_discons = len(discon_structs)

    def set_populated(self):
        self.populated = True

    def load_from_file(self):
        # check for initialization of file name
        if not self.file_name:
            raise Exception('Call to load_from_file() when file_name not initialized')
        
        # open the previously specified model file
        f_in = open(self.file_name, "r");
        
        # read in number of model parameters
        self.num_parameters = int(f_in.readline().strip())

        # read in volumetric parameter configurations (number of splines + descriptor)
        self.parameters = []
        for p in range(self.num_parameters):
            line = f_in.readline().strip().split()
            self.parameters.append(Parameter((int(line[0]), int(line[1]), line[2])))

        # read in number of peturbed discontinuities
        self.num_discons = int(f_in.readline().strip())

        # read in volumetric parameter configurations (number of splines + descriptor)
        self.discons = []
        for p in range(self.num_discons):
            line = f_in.readline().strip().split()
            self.discons.append(Discontinuity((int(line[0]), line[1])))

        # read in the vertical b-spline knot depths
        str = StringIO(f_in.readline())
        self.bspl_knots = np.loadtxt(str)

        # read in the remainder of the file (i.e. the model) as a numpy ndarray ...

        #  first for the volumetric parameters 
        for p in self.parameters:
            model_values = np.zeros((p.num_bspl, p.num_sspl))
            for l in range(p.num_bspl):
                strio = StringIO(f_in.readline())
                model_values[l,:] = np.loadtxt(strio)
            p.set_values(model_values)

        #  and the for the perturbed discontinuities
        for d in self.discons:
            strio = StringIO(f_in.readline())
            d.set_values(np.loadtxt(strio))

        # close the file
        f_in.close()

        # set the populated flag to True
        self.populated = True

    def save_to_file(self):
        # check for model population
        if not self.populated:
            raise Exception('Call to save_to_file() when model is not populated')

        # check for initialization of file name
        if not self.file_name:
            raise Exception('Call to save_to_file() when file_name not initialized')
        
        # open the previously specified model file
        f_out = open(self.file_name, "w");
        
        # write number of model parameters
        f_out.write('%i\n' % (self.num_parameters))

        # write in volumetric parameter configurations (number of splines + descriptor)
        for p in self.parameters:
            f_out.write('%i %i %s\n' % (p.num_sspl, p.num_bspl, p.desc))

        # write number of perturbed discontinuities
        f_out.write('%i\n' % (self.num_discons))

        # write in volumetric parameter configurations (number of splines + descriptor)
        for d in self.discons:
            f_out.write('%i %s\n' % (d.num_sspl, d.desc))
        
        # write out the b-spline knot radii
        np.savetxt(f_out, self.bspl_knots.reshape((1, self.bspl_knots.size)), fmt = '%f')

        # write out remainder of file using np.savetxt() ...

        # write out the parameter values
        for p in self.parameters:
            np.savetxt(f_out, p.values)

        # write out the discontinuity perturbations
        for d in self.discons:
            np.savetxt(f_out, d.values)

        # close the model file
        f_out.close()

