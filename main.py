#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
                                  .=""=.
                                |  d  b  |
                                \   /\   /
                               ,/'-=\/=-'\,
                              / /        \ \
                             | /          \ |
                             \/ \        / \/
                                 '.    .'
                                 _|`~~`|_
                                 / _  _ \
                                 /|\  /|\
              Nonadiabatic "on-the-fly" Molecular Dynamics.....
                       based on Zhu-Nakamura Theory
                              written by  python
"""
import copy
import os
import shutil
import sys
import time

import numpy as np

# All units in the code are a.u./hartree
# The units of initial coordination/velocity  is  Angstrom/Å  and Bohr / a.u.
np.set_printoptions(precision=8)
potential_energy = []  # molecular dynamics  energy
kin_energy = []  # molecular  kinetic energy
total_energy = []  # `kinetic + potential
q1 = []
q2 = []
q3 = []
coord = []  # x y z
momentum = []  # px py pz
rot_matrix = []
grad = []  #
scf_energy = []
hop_time = [0]  # hop time in dynamic programs
element = []
element_mass = []
atom_list = []
atom_value_range = []

# software 1.GAUSSIAN  2.ORCA  3.BDF 4.MOLPRO 5.MOLCAS
global soft
with open('run.inp', 'r') as f:
    for n, line in enumerate(f):
        if n == 7:
            soft = str(line.split()[0])
if soft == '1':
    from GAUSSIAN import *
elif soft == '2':
    from ORCA import *
elif soft == '3':
    from BDF import *
elif soft == '4':
    from MOLPRO import *
elif soft == '5':
    from MOLCAS import *
else:
    print('Soft is unknown')
    sys.exit()


def basis_constant():
    global eV, ang, fs, amu, pi, velo
    eV = 27.21138602
    ang = 0.529177257507  # Angstrom/Å e-10
    fs = 0.0241888439241  # 1 a.u. = 2.4188*10e-17
    amu = 1822.88853006  # Relative atomic mass
    pi = np.pi
    velo = 2.18769126364 * 10e6  # 1 a.u. = 2.187*10e4  m/s


def set_var():
    global SI_step_time, step_time, total_time, sim_begin_time, dynamics_time, initail_step_time
    global states_involved, dyn_states, states_max, states_min
    global nloop
    global soft
    global threshold
    global flag_gap,flag_rot 
    flag_gap = False
    flag_rot = False
    if os.path.exists('run.inp'):
        """
        states_involved : 4
        dyn_states : 3 -> (S3)
        states_max : 4 
        states_min : 1
        """
        with open('run.inp', 'r') as f:
            for n, line in enumerate(f):
                if n == 1:
                    line = [float(x) for x in line.split(",")]
                    sim_begin_time, SI_step_time, total_time = line[:3]
                    initail_step_time = SI_step_time
                    step_time = SI_step_time / fs
                elif n == 3:
                    nloop = int(line.split(",")[0])
                    dynamics_time = sim_begin_time
                elif n == 5:
                    line = [int(x) for x in line.split(",")]
                    states_involved, dyn_states, states_max, states_min = line[:4]
                elif n == 9:
                    threshold = float(line.split(",")[0])
                elif n == 11:
                    if line.split()[0] in ("Yy"):
                        flag_rot = True
    else:
        print("run.inp does not exist", flush=True)
        sys.exit(1)
    global inputfile, inputfile_d, inputfile_u
    global outputfile, outputfile_d, outputfile_u
    filename = None
    suffix = None
    suffix_1 = None
    if soft == '1':
        filename = 'gauss'
        suffix = '.gjf'
        suffix_1 = '.log'
    elif soft == '2':
        filename = 'orca'
        suffix = '.inp'
        suffix_1 = '.out'
    elif soft == '3':
        filename = 'bdf'
    elif soft == '4':
        filename = 'molpro'
    elif soft == '5':
        filename = 'molcas'
    else:
        pass
    inputfile = filename + suffix
    inputfile_d = filename + '_d' + suffix
    inputfile_u = filename + '_u' + suffix
    outputfile = filename + suffix_1
    outputfile_d = filename + '_d' + suffix_1
    outputfile_u = filename + '_u' + suffix_1


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


def current_time():
    return time.asctime(time.localtime(time.time()))


def get_initial_condition():
    def get_element_mass():
        global natom
        with open('initial_condition', 'r') as f:  # element x y z  p_x p_y p_z
            for value in f:
                if not value.isspace():
                    data = value.split()[0].capitalize()
                    element.append(data)
                    mass_au = float(masses[data]) * amu
                    element_mass.append(mass_au)
                else:
                    break
        natom = len(element)

    def get_key_element():
        if os.path.exists('geom.inp'):
            with open('geom.inp', 'r') as f:
                for line in f:
                    if not line.isspace():
                        if len(line.split()) == 8:
                            tmp = [int(i) for i in line.split()[
                                :4] if 0 < int(i) <= natom]
                            atom_list.append(tmp)
                            tmp_1 = [float(i) for i in line.split()[4:]]
                            atom_value_range.append(tmp_1)
                        else:
                            print(
                                "'geom.inp' format is error eg: atom1 atom2 atom3 atom4 num1  num2 num3 num4")
                            sys.exit(1)
                    else:
                        continue
        else:
            print("'geom.inp' does not exist")
            sys.exit(1)

    def get_position_momentum_matrix():
        # element x y z  p_x p_y p_z
        with open('initial_condition', 'r') as f:
            position_matrix = []
            momentum_matrix = []
            for value in f:
                if not value.isspace():
                    data = [float(i) / ang for i in value.split()[1:4]]
                    data1 = [float(i) for i in value.split()[-3:]]
                    position_matrix.append(data)
                    momentum_matrix.append(data1)
                else:
                    break
            return np.array(position_matrix), np.array(momentum_matrix)

    get_element_mass()
    get_key_element()
    xyz, mom = get_position_momentum_matrix()
    coord.append(xyz)
    momentum.append(mom)


def analyse_result(filename=None):
    """
    grad, potential_energy, total_energy
    """
    Filename = None
    if soft == '1':
        Filename = 'gauss.log'
    elif soft == '2':
        Filename = 'orca.out'
    elif soft == '3':
        Filename = 'bdf.out'
    elif soft == '4':
        Filename = 'molpro.out'
    elif soft == '5':
        Filename = 'molcas.out'
    else:
        pass
    if filename:
        Filename = filename
    grad.append(get_grad_matrix())
    E_exc, E_scf = get_energy()
    scf_energy.append(E_scf)
    if dyn_states >= 1 or states_involved <= 1:
        # S0 dynamics and the excited dynamics
        potential_energy.append(E_exc)
    else:  # after hopping S0 dynamics
        print("Begin calculate the excited states energy at %s" %
              current_time(), flush=True)
        read_wavefunction()
        # Get the keywords of the previous input file
        #  renew_calc_states(states_involved, inputfile) #error *****
        renew_calc_states(states_involved-1, inputfile,
                          remove="force")  # root=states_involved - 1
        software_running()
        E_exc, E_scf = get_energy()
        # get the excited stated excited energy
        potential_energy.append(E_exc)
        delete_wavefunction()
        #  renew_calc_states(dyn_states, inputfile) # error  *****
        renew_calc_states(dyn_states, inputfile, add="force")


def update_position_matrix(position_matrix, momentum_matrix, grad_matrix, element_mass):
    """
    x_n+1 = x_n +  p_n /m_i * Δt  + 1/2mi  * F_n * Δt^2
    """
    position_matrix_back = []
    for i in range(natom):
        data = position_matrix[i] + (momentum_matrix[i] / element_mass[i]) * \
            step_time - 0.5 * \
            (grad_matrix[i] / element_mass[i]) * step_time ** 2
        position_matrix_back.append(data)
    return np.array(position_matrix_back)


# notice atom force symbol positive + or negative  -


def update_momentum_matrix(momentum_matrix, grad_matrix, grad_matrix_back):
    """
    p_n+1 = p_n +　0.5 * (F_n+1 +F_n) * Δt
    """
    velocity_matrix_back = []
    for i in range(natom):
        data = momentum_matrix[i] - 0.5 * \
            (grad_matrix[i] + grad_matrix_back[i]) * step_time
        velocity_matrix_back.append(data)
    return np.array(velocity_matrix_back)


def critical_value(*args):
    def __dihedral_angle(p1, p2, p3, p4):
        q1 = np.subtract(p2, p1)
        q2 = np.subtract(p3, p2)
        q3 = np.subtract(p4, p3)
        # Calculate cross vectors
        q1_x_q2 = np.cross(q1, q2)
        q2_x_q3 = np.cross(q2, q3)
        # Calculate normal vectors
        n1 = q1_x_q2 / np.linalg.norm(q1_x_q2)
        n2 = q2_x_q3 / np.linalg.norm(q2_x_q3)
        # Calculate unit vectors
        u1 = n2
        u3 = q2 / np.linalg.norm(q2)
        u2 = np.cross(u3, u1)
        # Calculate cosine and sine
        cos_theta = np.dot(n1, u1)
        sin_theta = np.dot(n1, u2)
        theta = -np.arctan2(sin_theta, cos_theta)
        theta_deg = np.degrees(theta)
        return theta_deg

    def __bond_angle(p1, p2, p3):
        q1 = np.subtract(p2, p1)
        q2 = np.subtract(p3, p2)
        cos_theta = np.dot(q1, q2) / (np.linalg.norm(q1) * np.linalg.norm(q2))
        theta = np.arccos(cos_theta)
        theta_deg = 180.00 - np.degrees(theta)
        return theta_deg

    def __bond_length(p1, p2):
        length = np.linalg.norm(np.subtract(p2, p1))
        return length

    N = len(args)
    if N == 2:
        return __bond_length(*args)
    elif N == 3:
        return __bond_angle(*args)
    elif N == 4:
        return __dihedral_angle(*args)
    else:
        print('too many parameters')
        sys.exit()


def print_matrix(value):
    for i in range(natom):
        ele = format(element[i], '<5s')
        xyz = "".join(format(x, '>18.10f')for x in value[i])
        print(ele + xyz, flush=True)

def coord_mom_rot(coord1, coord2):
    """
    Kabsch algorithm
    url=https://en.wikipedia.org/wiki/Kabsch_algorithm
    """
    global flag_rot
    if not flag_rot: return np.eye(3)
    # natom = len(coord1)
    coord_1 = copy.deepcopy(coord1)
    coord_2 = copy.deepcopy(coord2)
    coord_1 -= np.mean(coord_1, axis=0)
    coord_2 -= np.mean(coord_2, axis=0)
    cov_matrix = coord_1.T @ coord_2
    U, sigma, VT = np.linalg.svd(cov_matrix)
    d = np.sign(np.linalg.det(U @ VT))
    diag = np.identity(3)
    diag[2][2] = d
    rot_mat= VT.T @ diag @ U.T
    # rot_coord = coord_2 @ rot_matrix
    # RMSD = np.sqrt(np.sum((coord_1 - rot_coord) ** 2) / natom)
    return  rot_mat 

def coord_coorections(geom, mass):
    """
    Modifie the coordinate to center of mass
    geom: N * 3;mass : 1 * N
    """
    mass_N = np.array(mass).reshape((len(mass), 1))
    coord_rot_matrix = coord_mom_rot(coord[0] * mass_N , geom * mass_N)
    new_geom = geom @ coord_rot_matrix
    return new_geom, coord_rot_matrix

def mom_coorections(mom1, mom2):
    """
    mom: N * 3
    """
    mom_rot_matrix = coord_mom_rot(mom1, mom2)
    new_mom = mom2 @ mom_rot_matrix
    return new_mom

def on_the_fly():  # 在software check——hooping 之后
    """
    renew the current step momentum and the nex step position
    Nuclear motion by velocity verlet Algorithm
    x_n+1 = x_n +  p_n /m_i * Δt  + 1/2mi * F_n * Δt^2
    p_n+1 = p_n +　0.5 * (F_n+1 +F_n)  * Δt
    """
    global q1, q2, q3
    if nloop == 0:  # first step
        # coord.append(update_position_matrix(coord[nloop], momentum[nloop], grad[nloop], element_mass))  
        # 2step position
        next_coord = update_position_matrix(coord[nloop], momentum[nloop], grad[nloop], element_mass)
        next_coord_c, rot = coord_coorections(next_coord, element_mass)
        coord.append(next_coord_c)
        rot_matrix.append(rot)
        q1 = coord[nloop]
        q2 = coord[nloop + 1]
    else:  # next step
        # the next position  & the current step momentum
        # momentum.append(update_momentum_matrix(momentum[nloop - 1], grad[nloop], grad[nloop - 1]))
        # coord.append(update_position_matrix(coord[nloop], momentum[nloop], grad[nloop], element_mass))
        next_mom = update_momentum_matrix(momentum[nloop - 1], grad[nloop], grad[nloop - 1])
        next_mom_c = next_mom @ rot_matrix[nloop-1]
        momentum.append(next_mom_c)
        next_coord = update_position_matrix(coord[nloop], momentum[nloop], grad[nloop], element_mass)
        next_coord_c, rot = coord_coorections(next_coord, element_mass)
        # next_mom_c = mom_coorections(momentum[0], next_mom)
        # print("This difference between next_mom and next_mom_c is %14.8f - %14.8f = %14.8f " %(calculate_kinetic(next_mom), calculate_kinetic(next_mom_c), calculate_kinetic(next_mom)-calculate_kinetic(next_mom_c)))
        coord.append(next_coord_c)
        rot_matrix.append(rot)
        q1 = coord[nloop - 1]
        q2 = coord[nloop]
        q3 = coord[nloop + 1]


def check_hopping():  # when nloop >=2,begin check hopping
    global dyn_states
    global nloop
    hop_type = None
    states = np.arange(states_involved)  # (0,1)
    # calculate energy difference of the q1 q2 q3 point　 suffix "u" is up  "d" is down

    if dyn_states == states[-1]:  # the highest state  e.g. S1 -> S0
        print("The current states(%d) is the highest state" %
              dyn_states, flush=True)
        state_d_num = dyn_states - 1
        delta_q1_u = 0.000000
        delta_q2_u = 0.000000
        delta_q3_u = 0.000000
        delta_q1_d = (potential_energy[nloop - 2][dyn_states] -
                      potential_energy[nloop - 2][state_d_num]) * eV
        delta_q2_d = (potential_energy[nloop - 1][dyn_states] -
                      potential_energy[nloop - 1][state_d_num]) * eV
        delta_q3_d = (potential_energy[nloop][dyn_states] -
                      potential_energy[nloop][state_d_num]) * eV
        print(
            "The energy differences between %s and %s states at the least three point(q1,q2,q3-current) are: %12.6f "
            "%12.6f %12.6f "
            % (dyn_states, state_d_num, delta_q1_d, delta_q2_d, delta_q3_d), flush=True)
    elif dyn_states == states[0]:  # the lowest state  e.g. S0 dynamic
        state_u_num = dyn_states + 1
        delta_q1_d = 0.00000
        delta_q2_d = 0.00000
        delta_q3_d = 0.00000
        if dyn_states == 0 and states_involved <= 1:  # S0 dynamic
            delta_q1_u = 0.00000
            delta_q2_u = 0.00000
            delta_q3_u = 0.00000
        else:
            delta_q1_u = (potential_energy[nloop - 2][state_u_num] -
                          potential_energy[nloop - 2][dyn_states]) * eV
            delta_q2_u = (potential_energy[nloop - 1][state_u_num] -
                          potential_energy[nloop - 1][dyn_states]) * eV
            delta_q3_u = (potential_energy[nloop][state_u_num] -
                          potential_energy[nloop][dyn_states]) * eV
        print("The current states(%d) is the lowest state" % dyn_states)
        if not (dyn_states == 0 and states_involved <= 1):  # S0 dynamic after the hopping
            print(
                "The energy differences between %s and %s states at the least three point(q1,q2,q3-current) are: %12.6f "
                "%12.6f %12.6f "
                % (dyn_states, state_u_num, delta_q1_u, delta_q2_u, delta_q3_u), flush=True)
    else:  # the middle state
        state_u_num = dyn_states + 1
        state_d_num = dyn_states - 1
        state_m_mum = dyn_states
        delta_q1_u = (potential_energy[nloop - 2][state_u_num] -
                      potential_energy[nloop - 2][state_m_mum]) * eV
        delta_q2_u = (potential_energy[nloop - 1][state_u_num] -
                      potential_energy[nloop - 1][state_m_mum]) * eV
        delta_q3_u = (potential_energy[nloop][state_u_num] -
                      potential_energy[nloop][state_m_mum]) * eV
        delta_q1_d = (potential_energy[nloop - 2][state_m_mum] -
                      potential_energy[nloop - 2][state_d_num]) * eV
        delta_q2_d = (potential_energy[nloop - 1][state_m_mum] -
                      potential_energy[nloop - 1][state_d_num]) * eV
        delta_q3_d = (potential_energy[nloop][state_m_mum] -
                      potential_energy[nloop][state_d_num]) * eV
        print("The current states(%d) is middle state" %
              dyn_states, flush=True)
        print(
            "The energy differences between %s and %s states at the least three point(q1,q2,q3-current) are: %12.6f "
            "%12.6f %12.6f "
            % (dyn_states, state_d_num, delta_q1_d, delta_q2_d, delta_q3_d), flush=True)
        print(
            "The energy differences between %s and %s states at the least three point(q1,q2,q3-current) are: %12.6f "
            "%12.6f %12.6f "
            % (dyn_states, state_u_num, delta_q1_u, delta_q2_u, delta_q3_u), flush=True)
    """
    check the type of hopping
    1. A double hop-upward and downward
    2. A downward hop
    3. A upward hop
    """
    if delta_q1_d > delta_q2_d and delta_q2_d < delta_q3_d:
        if delta_q1_u > delta_q2_u and delta_q2_u < delta_q3_u and delta_q2_u < kin_energy[nloop - 1]:
            if delta_q2_u <= threshold and delta_q2_d <= threshold:
                hop_type = 1
            elif delta_q2_d <= threshold < delta_q2_u:
                hop_type = 2
            elif delta_q2_u <= threshold < delta_q2_d:
                hop_type = 3
            else:
                print("No minimum energy separation smaller than %s ev was found at %d step " % (
                    threshold, nloop), flush=True)
        else:
            if delta_q2_d <= threshold:
                hop_type = 2
            else:
                print(
                    "The minimum energy separations for downward transitions are larger than  %s" % threshold,
                    flush=True)
    elif delta_q1_u > delta_q2_u and delta_q2_u < delta_q3_u and delta_q2_u < kin_energy[nloop - 1]:
        if delta_q2_u <= threshold:
            hop_type = 3
        else:
            print(" The minimum energy separations for upward transitions are larger \
                 than %s at %d step" % (threshold, nloop), flush=True)
    else:
        print("No local minimum between PESs observed at %d step" %
              nloop, flush=True)

    if hop_type == 2 and nloop >= (hop_time[-1] + 2):
        hop_direction = 'D'

        # calculate = q3 down_state grad
        # replace_coordinate(q3)  #the current coordinate is q3 ,so this is no significance
        renew_calc_states((dyn_states - 1), inputfile, inputfile_d)
        software_running(inputfile_d)
        grad_q3d = get_grad_matrix(outputfile_d)
        # calculate = q1  down_state grad

        replace_coordinate(q1)
        renew_calc_states((dyn_states - 1), inputfile, inputfile_d)
        software_running(inputfile_d)
        grad_q1d = get_grad_matrix(outputfile_d)

        hop_p, delta_grad_q2, mom_direction_factor, grad_q2d, grad_q2u = get_hop_factor(
            grad_q3d, grad[nloop], grad_q1d, grad[nloop - 2], hop_direction)
        rand_p = np.random.rand()
        print("The calculated hopping factor and random hopping factor is %14.6f %14.6f" % (
            hop_p, rand_p), flush=True)
        if hop_p >= rand_p:
            print("The downer Hopping succeed at %s" %
                  current_time(), flush=True)
            hopping_renew(delta_grad_q2, mom_direction_factor,
                          hop_direction, grad_q2d, grad_q2u)
        else:
            print("The downer Hopping failure %s" % current_time(), flush=True)
    elif hop_type == 3 and nloop >= (hop_time[-1] + 2):
        hop_direction = 'U'

        # calculate = q3 up_state grad
        renew_calc_states((dyn_states + 1), inputfile, inputfile_u)
        software_running(inputfile_u)
        grad_q3u = get_grad_matrix(outputfile_u)

        # calculate = q1  down_state grad
        replace_coordinate(q1)
        renew_calc_states((dyn_states + 1), inputfile, inputfile_u)
        software_running(inputfile_u)
        grad_q1u = get_grad_matrix(outputfile_u)

        hop_p, delta_grad_q2, mom_direction_factor, grad_q2d, grad_q2u = get_hop_factor(
            grad[nloop], grad_q3u, grad[nloop - 2], grad_q1u, hop_direction)
        rand_p = np.random.rand()
        print("The calculated hopping factor and random hopping factor is %14.6f %14.6f" % (
            hop_p, rand_p), flush=True)
        if hop_p >= rand_p:
            print("The upper Hopping succeed at %s" %
                  current_time(), flush=True)
            hopping_renew(delta_grad_q2, mom_direction_factor,
                          hop_direction, grad_q2d, grad_q2u)
        else:
            print("The upper Hopping failure at %s" %
                  current_time(), flush=True)
    elif hop_type == 1 and nloop >= (hop_time[-1] + 2):  # up median down
        grad_q3m = grad[nloop]
        grad_q1m = grad[nloop - 2]

        """
        calculate q3, q1 and q2 up state grad
        """
        hop_direction = 'U'
        hop_p = 0

        # calculate = q3 up_state grad
        renew_calc_states((dyn_states + 1), inputfile, inputfile_u)
        software_running(inputfile_u)
        grad_q3u = get_grad_matrix(outputfile_u)

        # calculate = q1 up_state grad
        replace_coordinate(q1)
        renew_calc_states((dyn_states + 1), inputfile, inputfile_u)
        software_running(inputfile_u)
        grad_q1u = get_grad_matrix(outputfile_u)

        hop_p_u, delta_grad_q2_u, mom_direction_factor_u, grad_q2m, grad_q2u = get_hop_factor(
            grad_q3m, grad_q3u, grad_q1m, grad_q1u, hop_direction)

        """
        calculate q3, q1 and q2 down state grad
        """
        hop_direction = 'D'
        # calculate = q3 down_state grad
        renew_calc_states((dyn_states - 1), inputfile, inputfile_d)
        software_running(inputfile_d)
        grad_q3d = get_grad_matrix(outputfile_d)

        # calculate = q1 down_state grad
        replace_coordinate(q1)
        renew_calc_states((dyn_states - 1), inputfile, inputfile_d)
        software_running(inputfile_d)
        grad_q1d = get_grad_matrix(outputfile_d)

        hop_p_d, delta_grad_q2_d, mom_direction_factor_d, grad_q2d, grad_q2m_1 = get_hop_factor(
            grad_q3d, grad_q3m, grad_q1d, grad_q1m, hop_direction)

        if hop_p_u > hop_p_d:
            hop_direction = 'U'
            rand_p = np.random.rand()
            print("The calculated upper/downer hopping factor and random hopping factor is %14.6f %14.6f %14.6f" %
                  (hop_p_u, hop_p_d, rand_p), flush=True)
            if hop_p_u >= rand_p:
                print("The upper Hopping succeed at %s" %
                      current_time(), flush=True)
                hopping_renew(delta_grad_q2_u, mom_direction_factor_u,
                              hop_direction, grad_q2m, grad_q2u)
            else:
                print("The upper Hopping failure at %s" %
                      current_time(), flush=True)
        else:
            hop_direction = 'D'
            rand_p = np.random.rand()
            print("The calculated upper/downer hopping factor and random hopping factor is %14.6f %14.6f %14.6f" %
                  (hop_p_u, hop_p_d, rand_p), flush=True)
            if hop_p_d >= rand_p:
                print("The downer Hopping succeed at %s" %
                      current_time(), flush=True)
                hopping_renew(delta_grad_q2_d, mom_direction_factor_d,
                              hop_direction, grad_q2d, grad_q2m_1)
            else:
                print("The downer Hopping failure %s" %
                      current_time(), flush=True)


def get_hop_factor(grad_q3d, grad_q3u, grad_q1d, grad_q1u, hop_direction):
    """
    Parameters:
    :param grad_q3d: the gradient of q3 down states
    :param grad_q3u: the gradient of q3 up states
    :param grad_q1d: the gradient of q1 down states
    :param grad_q1u: the gradient of q1 up states
    :param hop_direction: hopping direction ("D" or "U")
    :return: hop, delta_grad_q2, mom_direction_factor, grad_q2d, grad_q2u
    """

    grad_q2d = np.zeros((natom, 3))
    grad_q2u = np.zeros((natom, 3))
    mom_direction_factor = np.zeros((natom, 3))
    delta_grad_q2 = 0
    F12 = 0
    if hop_direction == 'D':  # downward  up
        Ex = (potential_energy[nloop][dyn_states] +
              potential_energy[nloop][dyn_states - 1]) * 0.50000
        Vx = (potential_energy[nloop][dyn_states] -
              potential_energy[nloop][dyn_states - 1]) * 0.50000
        for i in range(natom):
            # The coordinates of atom i at q1 and q3 have not changed
            for j in range(3):
                if np.abs((q3[i][j] - q1[i][j])) <= 0.00005:
                    grad_q2u[i][j] = grad_q2d[i][j] = 0
                else:
                    # The current states is upper states =>q2u
                    grad_q2d[i][j] = - (1 / (q3[i][j] - q1[i][j])) * (grad_q3d[i][j] * (q2[i][j] - q1[i][j])
                                                                      - grad_q1u[i][j] * (q2[i][j] - q3[i][j]))
                    grad_q2u[i][j] = - (1 / (q3[i][j] - q1[i][j])) * (grad_q3u[i][j] * (q2[i][j] - q1[i][j])
                                                                      - grad_q1d[i][j] * (q2[i][j] - q3[i][j]))
    elif hop_direction == 'U':
        Ex = (potential_energy[nloop][dyn_states] +
              potential_energy[nloop][dyn_states + 1]) * 0.50000
        Vx = (potential_energy[nloop][dyn_states] -
              potential_energy[nloop][dyn_states + 1]) * 0.50000
        for i in range(natom):
            # Coordinate difference between q1 and q3
            for j in range(3):
                if np.abs((q3[i][j] - q1[i][j])) <= 0.00005:
                    grad_q2u[i][j] = grad_q2d[i][j] = 0
            else:
                # The current states is downer states =>q2d
                grad_q2d[i][j] = -(1 / (q3[i][j] - q1[i][j])) * (grad_q3u[i][j] * (q2[i][j] - q1[i][j])
                                                                 - grad_q1d[i][j] * (q2[i][j] - q3[i][j]))
                grad_q2u[i][j] = -(1 / (q3[i][j] - q1[i][j])) * (grad_q3d[i][j] * (q2[i][j] - q1[i][j])
                                                                 - grad_q1u[i][j] * (q2[i][j] - q3[i][j]))

    for i in range(natom):
        for j in range(3):
            # \sum{\frac{(F^2_i - F^1_i)^2}{m_i}}
            delta_grad_q2 += (grad_q2u[i][j] -
                              grad_q2d[i][j]) ** 2 / element_mass[i]
            # \frac{F^2_i - F^1_1 }{\sqrt{m_i}}
            mom_direction_factor[i][j] = (
                grad_q2u[i][j] - grad_q2d[i][j]) / np.sqrt(element_mass[i])
            # \sum{F^1_i * F^2_i}
            F12 += grad_q2d[i][j] * grad_q2u[i][j]

    # calculate a^2 and b2
    f_aa = delta_grad_q2 / (16 * Vx ** 3)  # a.u. Reduced Planck constant =1
    f_bb = (total_energy[nloop - 1] - Ex) / (2 * Vx)
    if f_aa > 1000:
        hop_p = 1.0
    elif f_aa < 0.001:
        hop_p = 0.0
    else:
        if F12 >= 0:
            hop_p = np.exp(-np.pi / (4 * np.sqrt(f_aa)) *
                           np.sqrt(2 / (f_bb + np.sqrt(f_bb ** 2 + 1))))
        else:
            hop_p = np.exp(-np.pi / (4 * np.sqrt(f_aa)) *
                           np.sqrt(2 / (f_bb + np.sqrt(f_bb ** 2 - 1))))

    # print key value
    print("**************")
    if hop_direction == "D":
        print("The current(upper) states %s " % dyn_states, flush=True)
        print("The downer states %s" % (dyn_states-1), flush=True)
    if hop_direction == "U":
        print("The current(downer) states %s" %
              dyn_states, flush=True)
        print("The upper states is states %s" % (dyn_states+1), flush=True)
    print("Coordinates of q1 (a.u.)")
    print_matrix(q1)
    print("Gradient of q1 in downer states ", flush=True)
    print_matrix(grad_q1d)
    print("Gradient of q1 in upper states", flush=True)
    print_matrix(grad_q1u)
    print("Coordinates of q3 (a.u.)")
    print_matrix(q3)
    print("Gradient of q3 in downer states", flush=True)
    print_matrix(grad_q3d)
    print("Gradient of q3 in upper states", flush=True)
    print_matrix(grad_q3d)
    print("Coordinates of q2 (a.u.)")
    print_matrix(q2)
    print("Gradient of q2 in downer states")
    print_matrix(grad_q2d)
    print("Gradient of q2 in upper states")
    print_matrix(grad_q2u)
    print("**************")
    print("The (F2-F1)^2/m ,F12 is %18.10f %18.10f" %
          (delta_grad_q2, F12), flush=True)
    print("The total-energy, Vx, Ex is %18.10f %18.10f %18.10f" %
          (total_energy[nloop-1], Vx, Ex), flush=True)
    print("a^2, b^2, hop_p is %18.10f %18.10f %18.10f" %
          (f_aa, f_bb, hop_p), flush=True)

    return hop_p, delta_grad_q2, mom_direction_factor, grad_q2d, grad_q2u


def hopping_renew(delta_grad_q2, mom_direction_factor, hop_direction, grad_q2d, grad_q2u):
    """
    After the hopping, revise q2 momentum,kinetic energy and renew q3 coordinate
    """
    # Calculate every atom normalized momentum factor
    # s_i = \frac{\frac{F^2_i(q^2)-F^1_i(q^2)}{\sqrt{m_i}}}{\sqrt{\sum{\frac{(F^2_i-F^1_i)^2}{m_i}}}}
    # hop direction
    # n_i = \frac{s_i}{|s_i|}
    # increment factor k
    # k = \sqrt{1+\frac{U_+(q^2)-U_-(q^2)}{\sum_{i=1}^{N}\frac{P_{i//}^2(+)}{2m_i}}}

    global dyn_states, nloop
    hop_direction_normal = np.zeros((natom, 3))
    hop_direction_matrix = mom_direction_factor / np.sqrt(delta_grad_q2)  # s_i
    kinetic_parallel = 0
    if hop_direction == 'D':
        last_states = dyn_states - 1
    elif hop_direction == 'U':
        last_states = dyn_states + 1
    else:
        print("The hop direction is error")
        sys.exit()

    # delta_energy_q2 = potential_energy_q2(+) - potential_energy_q2(-)
    delta_energy_q2 = potential_energy[nloop - 1][dyn_states] - \
        potential_energy[nloop - 1][last_states]
    mom_parallel = np.zeros((natom, 3))
    mom_renew = np.zeros((natom, 3))

    # before the hopping
    for i in range(natom):
        # n_i Normalized Si
        hop_direction_normal[i] = hop_direction_matrix[i] / \
            np.linalg.norm(hop_direction_matrix[i])
        # p_parallel = (n_i . p_i) * n_i
        mom_parallel[i] = np.dot(
            hop_direction_normal[i], momentum[nloop - 1][i]) * hop_direction_normal[i]
        # Parallel kinetic energy
        kinetic_parallel += 0.50 * \
            np.linalg.norm(mom_parallel[i]) ** 2 / element_mass[i]

    if hop_direction == 'D':
        # higher state : e.g. S1->S0 k > 1
        increment_k = np.sqrt(1 + delta_energy_q2 / kinetic_parallel)
    else:
        # hop_direction == 'U': Lower state - higher state e.g. :S1->S2 k < 1
        increment_k = np.sqrt(1 - delta_energy_q2 / kinetic_parallel)
    # calculate the momentum after the hopping
    for i in range(natom):
        mom_renew[i] = momentum[nloop - 1][i] + (increment_k - 1) * \
            np.dot(momentum[nloop - 1][i],
                   hop_direction_normal[i]) * hop_direction_normal[i]
    print("The increment factor(k) for momentum adjustment is %.8f" %
          increment_k, flush=True)
    print("***************")
    # print("Normalized momentum direction factor(S_i) is ", flush=True)
    # print_matrix(mom_direction_factor)
    print("Normalized momentum direction factor(S_i) ", flush=True)
    print_matrix(hop_direction_matrix)
    print("Renormalized momentum direction factor (n_i)", flush=True)
    print_matrix(hop_direction_normal)
    print("The q2 old momentum is", flush=True)
    print_matrix(momentum[nloop - 1])
    print("The parallel momentum is %.8f eV" %
          (kinetic_parallel*eV), flush=True)
    print_matrix(mom_parallel)
    print("The q2 new momentum is ", flush=True)
    print_matrix(mom_renew)
    print("***************")

    before_dyn_states = dyn_states
    if hop_direction == 'D':
        dyn_states = dyn_states - 1
    elif hop_direction == 'U':
        dyn_states = dyn_states + 1

    # revise the dyn_states
    renew_calc_states(dyn_states, inputfile)

    # delete q2, q3 grad
    grad.pop()
    grad.pop()
    new_grad_q2 = np.zeros((natom, 3))
    # renew q2 grad
    if hop_direction == 'D':
        new_grad_q2 = grad_q2d
    elif hop_direction == 'U':
        new_grad_q2 = grad_q2u
    grad.append(new_grad_q2)

    # check  q2 kinetic energy before and after hopping
    q2_old_k = calculate_kinetic(momentum[nloop - 1])
    # delete q2 momentum and renew q2 momentum
    momentum.pop()
    mom_renew_c = mom_renew @ rot_matrix[nloop-1] #??????
    momentum.append(mom_renew_c)
    q2_new_k = calculate_kinetic(momentum[nloop - 1])

    print("The parallel kinetic is %.8f eV" %
          (kinetic_parallel*eV), flush=True)
    print("The q2 old/new kinetic energy is %15.10f eV %15.10f eV" %
          (q2_old_k*eV, q2_new_k*eV), flush=True)
    print("The kinetic energy difference is %15.8f eV" %
          ((q2_old_k - q2_new_k)*eV), flush=True)
    print("The energy difference is %15.10f eV" %
          (delta_energy_q2*eV), flush=True)
    total_old = potential_energy[nloop - 1][before_dyn_states] + q2_old_k
    total_new = potential_energy[nloop - 1][dyn_states] + q2_new_k
    print("The q2 old/new kinetic energy is %15.10f %15.10f " %
          (q2_old_k, q2_new_k), flush=True)
    print("The q2 old/new scf energy is %15.10f %15.10f" %
          (potential_energy[nloop - 1][before_dyn_states], potential_energy[nloop - 1][dyn_states]), flush=True)
    print("The q2 old/new total energy is %15.10f %15.10f" %
          (total_old, total_new), flush=True)

    # delete q3 coord
    q3_old = coord.pop()
    # delete q2 kinetic_energy /total_energy
    kin_energy.pop()
    kin_energy.append(q2_new_k)
    # delete q3 potential
    potential_energy.pop()
    # total_energy is constant during the hopping

    # delete q3 scf_energy/
    scf_energy.pop()  # ??????

    # deleta q3 rot_matrix 
    rot_matrix.pop() 

    # revise q2 kinetic-energy potential-energy total energy
    # delete q2 coordinate and momentum
    with open('simulation.xyz', 'r') as f,  open("%s.bak" % 'simulation.xyz', 'w+') as f_new:
        lines = f.readlines()
        f_new.writelines(lines[:-(natom + 2)])
        f_new.write(str(natom) + '\n')
        f_new.write('simulation time: t  =' +
                    format(dynamics_time - SI_step_time, '>10.2f') + ' \n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            xyz = ''.join(format(x * ang, '>18.10f')
                          for x in coord[nloop - 1][i])
            mom = ''.join(format(x, '>18.10f')
                          for x in momentum[nloop - 1][i])
            f_new.write(ele + xyz + mom + '\n')
    os.remove('simulation.xyz')
    os.rename("%s.bak" % 'simulation.xyz', 'simulation.xyz')

    # record hop time (q2, nloop -1)
    with open('traj_hopping.log', 'a+') as f:
        f.write('{:<8.2f}fs{:>5s} => {:<5s}'.format(dynamics_time - SI_step_time,
                                                    str(before_dyn_states), str(dyn_states)) + '\n')
        f.write('{:<15.8f}{:>15.8f}{:>15.8f}{:>15.8f}'.format(total_energy[nloop - 1],
                                                              potential_energy[nloop -
                                                                               1][before_dyn_states],
                                                              potential_energy[nloop -
                                                                               1][dyn_states],
                                                              kin_energy[nloop - 1]))
        f.write('   ' + '\n')

    # revise q2 traj_energy.log (q2) q3 is current coordinate
    with open('traj_energy.log', 'r') as f, open("%s.bak" % 'traj_energy.log', 'w+') as f_new:
        lines = f.readlines()
        f_new.writelines(lines[:-1])
        energy_involved = ''.join(
            format(potential_energy[nloop - 1][i], '>18.8f') for i in range(states_involved))
        test = '{:<8.2f}{:>18.8f}{:>18.8f}{:>18.8f}'.format(dynamics_time - SI_step_time, total_energy[nloop - 1],
                                                            potential_energy[nloop - 1][dyn_states], kin_energy[nloop - 1])
        f_new.write(test + '     ' + energy_involved + '\n')
    os.remove('traj_energy.log')
    os.rename("%s.bak" % 'traj_energy.log', 'traj_energy.log')

    # recalculate q3 gradient potential_energy gradient q2 coordinate, gradient momentum
    # coord.append(update_position_matrix(
    #     coord[nloop - 1], momentum[nloop - 1], grad[nloop - 1], element_mass))
    # coord.append(update_position_matrix(q2, mom_renew, new_grad_q2, element_mass))
    q3 = update_position_matrix(q2, mom_renew, new_grad_q2, element_mass) 
    q3_c, rot = coord_coorections(q3, element_mass)
    coord.append(q3_c)
    rot_matrix.append(rot)
    replace_coordinate(coord[nloop])
    software_running()
    analyse_result()

    # q3_kin = calculate_kinetic(update_momentum_matrix(
    #   momentum[nloop - 1], grad[nloop], new_grad_q2))
    q3_kin = calculate_kinetic(update_momentum_matrix(
        mom_renew, grad[nloop], new_grad_q2))
    q3_scf = scf_energy[nloop]
    q3_total = q3_kin + q3_scf

    print("The q3 kinetic/scf energy is %18.8f eV %18.10f" %
          (q3_kin*eV, q3_scf), flush=True)
    print("The q3 total energy is %18.10f" % (q3_total), flush=True)
    print("*************", flush=True)
    print("The q3 old coordinate(a.u.) is ", flush=True)
    print_matrix(q3_old)
    print("The q3 adjusting coordinate(a.u.) is", flush=True)
    print_matrix(coord[nloop])
    print("The q3 grad is ",flush=True)
    print_matrix(grad[nloop])
    print("*************", flush=True)


def check_result():
    # calculate bond length angle  dihedral
    global SI_step_time, step_time, flag_gap
    bond1 = [25, 28]
    bond2 = [22, 28]
    bond3 = [22, 31]
    bond4 = [31, 30]
    key1 = [coord[nloop][i - 1] * ang for i in bond1]
    key2 = [coord[nloop][i - 1] * ang for i in bond2]
    key3 = [coord[nloop][i - 1] * ang for i in bond3]
    key4 = [coord[nloop][i - 1] * ang for i in bond4]
    OH1 = critical_value(*key1)
    OH2 = critical_value(*key2)
    OH3 = critical_value(*key3)
    OH4 = critical_value(*key4)
    if not flag_gap:
        if np.abs(OH1 - OH2) <= 0.100 or np.abs(OH3 - OH4) <= 0.100:
            SI_step_time = 0.2
            step_time = SI_step_time / fs
            flag_gap = True
    else:
        if np.abs(OH1 - OH2) > 0.100 and np.abs(OH3 - OH4) > 0.100:
            SI_step_time = initail_step_time
            step_time = SI_step_time / fs
            flag_gap = False


def MO_draw():
    if dynamics_time <= 50.001 and nloop % 2 == 0:
        filename1 = '62_' + format(dynamics_time, '.1f') + '.cub'
        filename2 = '63_' + format(dynamics_time, '.1f') + '.cub'
        os.system("formchk gauss.chk &> /dev/null ")
        os.system("cubegen 0 mo=62 gauss.fchk %s 0 h &> /dev/null " %
                  filename1)
        os.system("cubegen 0 mo=63 gauss.fchk %s 0 h &> /dev/null " %
                  filename2)
        if not os.path.isdir('MO'):
            os.mkdir('MO')
        shutil.move(filename1, './MO')
        shutil.move(filename2, './MO')
        # os.system("ls MO &> /dev/null || mkdir MO && mv  %s %s MO/ " % (filename1, filename2))


def calculate_keyvalue():
    value = []
    for label in atom_list:
        c = [coord[nloop][i] * ang for i in label]
        a = "{:.8f}".format(critical_value(*c))
        value.append(a)
    return value


def calculate_kinetic(mom: list):
    """
    :param mom: the momentum of
    :return: the kinetic energy
    """
    kine = 0.0
    for i in range(natom):
        kine += 0.5 * np.sum(mom[i] ** 2) / element_mass[i]
    return kine


def print_result():
    # print simulation.xyz
    # element  x y z px py pz .10f
    E_kine = calculate_kinetic(momentum[nloop])
    kin_energy.append(E_kine)
    total_energy.append(E_kine + scf_energy[nloop])
    print("The total/kinetic/potential energy(a.u.) %16.8f %16.8f %16.8f at the %6.2f fs"
          % (total_energy[nloop], kin_energy[nloop], scf_energy[nloop], dynamics_time), flush=True)
    with open('simulation.xyz', 'a+') as f:
        f.write(str(natom) + '\n')
        f.write('simulation time: t  =' +
                format(dynamics_time, '>10.2f') + ' \n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            xyz = ''.join(format(x * ang, '>18.10f')
                          for x in coord[nloop][i])
            mom = ''.join(format(x, '>18.10f')
                          for x in momentum[nloop][i])
            f.write(ele + xyz + mom + '\n')

    # print traj_energy.log
    # simulation_time total_energy  potential_energy kinetic_energy  energy_involved  .8f
    with open('traj_energy.log', 'a+') as f:
        energy_involved = ''.join(
            format(potential_energy[nloop][i], '>18.8f') for i in range(states_involved))
        test = '{:<8.2f}{:>18.8f}{:>18.8f}{:>18.8f}'.format(dynamics_time, total_energy[nloop],
                                                            potential_energy[nloop][dyn_states], kin_energy[nloop])
        f.write(test + '     ' + energy_involved + '\n')

    # print traj_cicoe.log
    # print key coordinate  angle bond dihedral
    print_traj_cicoe(dynamics_time)
    with open('traj_coord.log', 'a+') as f:
        f.write('{:<12.2f}'.format(dynamics_time))
        for label in atom_list:
            c = [coord[nloop][i - 1] * ang for i in label]
            f.write('{:>12.4f}'.format(critical_value(*c)))
        f.write('\n')


def main():
    print("The Nonadiabatic on-the-fly Molecular Dynamics based on Zhu-Nakamura Theory", flush=True)
    print("--------------------------------------\n", flush=True)
    global dynamics_time, nloop
    basis_constant()
    set_var()
    check_initial(dyn_states)
    get_initial_condition()  # read coordinates from initial condition
    replace_coordinate(coord[0])
    while dynamics_time <= total_time:
        print("Trajectory calculation of step %s loop start at the %5.2f fs (%s)" % (
            nloop, dynamics_time, current_time()), flush=True)
        software_running()
        analyse_result()
        if nloop >= 2:
            check_hopping()
        # check_result()
        on_the_fly()
        print_result()
        # MO_draw()
        replace_coordinate(coord[nloop + 1])
        print("Trajectory calculation of step %s loop end at the %5.2f fs (%s)\n\n" % (
            nloop, dynamics_time, current_time()), flush=True)
        nloop += 1
        dynamics_time += SI_step_time
    print("The dynamics program has ended at the %5.2f fs (%s)\n\n" %
          (dynamics_time, current_time()), flush=True)


if __name__ == "__main__":
    main()

