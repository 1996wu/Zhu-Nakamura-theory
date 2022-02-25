import configparser
import os
import sys
import numpy as np

ang = 0.52917725750
amu = 1822.88853006  # Relative atomic mass
fs = 0.0241888439241

# Relative atomic mass
masses = {'X': 0, 'Ac': 227.028, 'Al': 26.981539, 'Am': 243, 'Sb': 121.757, 'Ar': 39.948, 'As': 74.92159, 'At': 210,
          'Ba': 137.327, 'Bk': 247, 'Be': 9.012182, 'Bi': 208.98037, 'Bh': 262, 'B': 10.811, 'Br': 79.904,
          'Cd': 112.411, 'Ca': 40.078, 'Cf': 251, 'C': 12.011, 'Ce': 140.115, 'Cs': 132.90543, 'Cl': 35.4527,
          'Cr': 51.9961, 'Co': 58.9332, 'Cu': 63.546, 'Cm': 247, 'Db': 262, 'Dy': 162.5, 'Es': 252, 'Er': 167.26,
          'Eu': 151.965, 'Fm': 257, 'F': 18.9984032, 'Fr': 223, 'Gd': 157.25, 'Ga': 69.723, 'Ge': 72.61,
          'Au': 196.96654, 'Hf': 178.49, 'Hs': 265, 'He': 4.002602, 'Ho': 164.93032, 'H': 1.00794, 'In': 114.82,
          'I': 126.90447, 'Ir': 192.22, 'Fe': 55.847, 'Kr': 83.8, 'La': 138.9055, 'Lr': 262, 'Pb': 207.2, 'Li': 6.941,
          'Lu': 174.967, 'Mg': 24.305, 'Mn': 54.93805,
          'Mt': 266, 'Md': 258, 'Hg': 200.59, 'Mo': 95.94, 'Nd': 144.24, 'Ne': 20.1797, 'Np': 237.048, 'Ni': 58.6934,
          'Nb': 92.90638, 'N': 14.00674, 'No': 259, 'Os': 190.2, 'O': 15.9994, 'Pd': 106.42, 'P': 30.973762,
          'Pt': 195.08, 'Pu': 244, 'Po': 209, 'K': 39.0983, 'Pr': 140.90765, 'Pm': 145, 'Pa': 231.0359, 'Ra': 226.025,
          'Rn': 222, 'Re': 186.207, 'Rh': 102.9055, 'Rb': 85.4678, 'Ru': 101.07, 'Rf': 261, 'Sm': 150.36,
          'Sc': 44.95591, 'Sg': 263,
          'Se': 78.96, 'Si': 28.0855, 'Ag': 107.8682, 'Na': 22.989768, 'Sr': 87.62, 'S': 32.066, 'Ta': 180.9479,
          'Tc': 98, 'Te': 127.6, 'Tb': 158.92534, 'Tl': 204.3833, 'Th': 232.0381, 'Tm': 168.93421, 'Sn': 118.71,
          'Ti': 47.88, 'W': 183.85, 'U': 238.0289, 'V': 50.9415, 'Xe': 131.29, 'Yb': 173.04, 'Y': 88.90585, 'Zn': 65.39,
          'Zr': 91.224}

config = configparser.ConfigParser()
config.read("run.ini")

sim_begin_time = config.getfloat('Time', 'Initialize_time')
SI_step_time = config.getfloat('Time', 'Delta_time')
total_time = config.getfloat('Time', 'Total_time')
dynamics_time = sim_begin_time
nloop = config.getint('Time', "Current_loop")


states_involved = config.getint('State', 'Total_state')
dyn_states = config.getint('State', 'Current_state')
states_max = config.getint('State', 'Max_state')
states_min = config.getint('State', 'Min_state')


soft = config.get('Run_soft', 'Soft').capitalize()
try:
    soft_path = config.get('Soft_path', soft + '_path')
except NoOptionError:
    pass

threshold = config.getfloat('Delta_energy', 'threshold')

flag_rot = config.getboolean('Const_cent', 'flag_rot')

path_initial = config.get("Initial_condition", 'coord')
natom = config.getint("Initial_condition", 'natom')

path_geom = config.get("Geom", 'atom_list')


def get_position_momentum_matrix(filename):
    # element x y z  p_x p_y p_z
    with open(filename, 'r') as f:
        position_matrix = []
        momentum_matrix = []
        for value in f:
            if not value.isspace():
                data = [float(i) / ang for i in value.split()[1:4]]
                data1 = [float(i) for i in value.split()[-3:]]
                position_matrix.append(data)
                momentum_matrix.append(data1)
        return np.array(position_matrix), np.array(momentum_matrix)


initial_coord, initial_mom = get_position_momentum_matrix(path_initial)


def get_element_mass(filename):
    em = []
    e = []
    with open(filename, 'r') as f:  # element x y z  p_x p_y p_z
        for value in f:
            if not value.isspace():
                data = value.split()[0].capitalize()
                e.append(data)
                mass_au = float(masses[data]) * amu
                em.append(mass_au)
    return len(e), np.array(em), e


natom, element_mass, element = get_element_mass(path_initial)


def get_key_element(filename):
    al = []
    avr = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if not line.isspace():
                    if len(line.split()) == 8:
                        tmp = [int(i) for i in line.split()[
                                               :4] if 0 < int(i) <= natom]
                        al.append(tmp)
                        tmp_1 = [float(i) for i in line.split()[4:]]
                        avr.append(tmp_1)
                    else:
                        raise Exception(
                            "'geom.inp' format is error eg: atom1 atom2 atom3 atom4 num1  num2 num3 num4")
    except FileNotFoundError:
        raise Exception("'geom.inp' does not exist")
    return al, avr


atom_list, atom_value_range = get_key_element(path_geom)

input_prefix = {"Gaussian": "gauss", "Orca": "orca", "Molpro": "molpro"}
input_suffix = {"Gaussian": "gjf", "Orca": "inp", "Molpro": "in"}
output_suffix = {"Gaussian": "log", "Orca": "out", "Molpro": "out"}
wfn_suffix = {"Gaussian": "chk", "Orca": "gbw", "Molpro": "wfu"}


def filename_var():
    filename = input_prefix[soft]
    suffix = input_suffix[soft]
    suffix_1 = output_suffix[soft]
    suffix_2 = wfn_suffix[soft]
    input_m = filename + "." + suffix
    input_d = filename + '_d' + "." + suffix
    input_u = filename + '_u' + "." + suffix
    output = filename + "." + suffix_1
    output_d = filename + '_d' + "." + suffix_1
    output_u = filename + '_u' + "." + suffix_1
    wfn = filename + "." + suffix_2
    return input_m, input_d, input_u, output, output_d, output_u, wfn


inputfile, inputfile_d, inputfile_u, outputfile, outputfile_d, outputfile_u, wfnfile = filename_var()


class EnterDir:
    def __init__(self, c):
        self.dir = c

    def __enter__(self):
        try:
            os.chdir(self.dir)
        except FileNotFoundError:
            sys.exit()

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir("../")

def hop_record():
    try:
        os.mkdir("HoppingRecord")
    except FileExistsError:
        pass 
hop_record()