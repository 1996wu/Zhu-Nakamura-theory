#!/home/wzb/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
                                  .=""=.
                                 / _  _ \
                                |  d  b  |
                                \   /\   /
                               ,/'-=\/=-'\,
                              / /        \ \
                             | /          \ |
                             \/ \        / \/
                                 '.    .'
                                 _|`~~`|_
                                 /|\  /|\
              Nonadiabatic "on-the-fly" Molecular Dynamics.....
                       based on Zhu-Nakamura Theory
                              written by  python
"""

import os
import re
import numpy as np
import sys
import time

from numpy.lib.function_base import diff 

# All units in the code are a.u./hartree
# The units of initial coordination/velocity  is ang/Å

eV = 27.21138602
ang = 0.529177257507  # Angstrom/Å e-10
fs = 0.0241888439241  # 1 a.u. = 2.4188*10e-17 s
amu = 1822.88853006  # Relative atomic mass
pi = 3.14159265358979323846
Velo = 2.18769126364*10e6  # 1 a.u. = 2.187*10e4  m/s

SI_step_time = 0.50  # unit:fs e-14
step_time = SI_step_time / fs
natom = 26
total_time = 1000 #unit fs
nloop = 0
threshold=0.3   ###energy difference between two states (unit eV)


dynamic_energy = [ ]   ## molecular dynamics  energy

print(Velo)
# Relative atomic mass
masses = {'X': 0, 'Ac': 227.028, 'Al': 26.981539, 'Am': 243, 'Sb': 121.757, 'Ar': 39.948, 'As': 74.92159, 'At': 210, 'Ba': 137.327, 'Bk': 247, 'Be': 9.012182, 'Bi': 208.98037, 'Bh': 262, 'B': 10.811, 'Br': 79.904, 'Cd': 112.411, 'Ca': 40.078, 'Cf': 251, 'C': 12.011, 'Ce': 140.115, 'Cs': 132.90543, 'Cl': 35.4527, 'Cr': 51.9961, 'Co': 58.9332, 'Cu': 63.546, 'Cm': 247, 'Db': 262, 'Dy': 162.5, 'Es': 252, 'Er': 167.26,
          'Eu': 151.965, 'Fm': 257, 'F': 18.9984032, 'Fr': 223, 'Gd': 157.25, 'Ga': 69.723, 'Ge': 72.61, 'Au': 196.96654, 'Hf': 178.49, 'Hs': 265, 'He': 4.002602, 'Ho': 164.93032, 'H': 1.00794, 'In': 114.82, 'I': 126.90447, 'Ir': 192.22, 'Fe': 55.847, 'Kr': 83.8, 'La': 138.9055, 'Lr': 262, 'Pb': 207.2, 'Li': 6.941, 'Lu': 174.967, 'Mg': 24.305, 'Mn': 54.93805,
          'Mt': 266, 'Md': 258, 'Hg': 200.59, 'Mo': 95.94, 'Nd': 144.24, 'Ne': 20.1797, 'Np': 237.048, 'Ni': 58.6934, 'Nb': 92.90638, 'N': 14.00674, 'No': 259, 'Os': 190.2, 'O': 15.9994, 'Pd': 106.42, 'P': 30.973762, 'Pt': 195.08, 'Pu': 244, 'Po': 209, 'K': 39.0983, 'Pr': 140.90765, 'Pm': 145, 'Pa': 231.0359, 'Ra': 226.025, 'Rn': 222, 'Re': 186.207, 'Rh': 102.9055, 'Rb': 85.4678, 'Ru': 101.07, 'Rf': 261, 'Sm': 150.36, 'Sc': 44.95591, 'Sg': 263,
          'Se': 78.96, 'Si': 28.0855, 'Ag': 107.8682, 'Na': 22.989768, 'Sr': 87.62, 'S': 32.066, 'Ta': 180.9479, 'Tc': 98, 'Te': 127.6, 'Tb': 158.92534, 'Tl': 204.3833, 'Th': 232.0381, 'Tm': 168.93421, 'Sn': 118.71, 'Ti': 47.88, 'W': 183.85, 'U': 238.0289, 'V': 50.9415, 'Xe': 131.29, 'Yb': 173.04, 'Y': 88.90585, 'Zn': 65.39, 'Zr': 91.224}





def software_running():
    current_time =  time.asctime( time.localtime(time.time()) )
    if os.path.isfile ("gaussian.sh"):
        print(" The %d time Gaussian begin running  at %s \n" %(nloop, current_time))
        os.system(" sh gaussian.sh ")
    else:
        print( "Gaussian  script is not found \n dynamic program has end at %s"  %current_time)
        sys.exit()

def  check_initial():
    print(1)

def replace_coordinate(filename,filename_new,new_coordinate):
    regex  = re.compile('[A-Za-z]{1,2}\s*(\s*(-?[0-9]+\.[0-9]*)){3}')
    regex_1 = re.compile('(\s*(-?[0-9]+\.[0-9]*)){3}')
    count = 0
    with open (filename, 'r') as f:
        with open (filename_new, 'w+') as f_new:
            for n , line in enumerate (f):
                if not regex.findall(line):
                    f_new.write(line)
                else:
                    replace_coord = '\t' + '\t'.join(str(i) for i in new_coordinate[count][:3])
                    f_new.write(re.sub(regex_1,replace_coord,line))
                    count +=1


def get_grad_matrix(filename,natom):
    regex = re.compile('Forces \(Hartrees/Bohr\)')
    data = [ ]
    addflag = False
    count = 0
    with open(filename,'r') as f:
        for line in f:
            if not regex.search(line) or addflag:
                    if addflag:
                        if count < 2:
                            count += 1 
                        else:
                            data.append(list(map(float,line.split()[2:])))
                            if (len(data)  == natom ): ##mathch line 
                                break
            else:
                addflag = True
    gradient_matrix = -np.array(data) #the gradient is opposite to the direction of force
    return gradient_matrix

def get_energy(filename):
    regex = re.compile('SCF Done')
    regex_1 = re.compile('Excited State\s*[0-9]*:')
    energy = [ ]
    E_ground = None
    with open (filename,'r') as f:
        for line in f:
            if regex.search(line):
                E_ground = float(line.split()[4])
                energy.append(E_ground)
            if regex_1.search(line): 
                energy.append(float(line.split()[4]) / eV + E_ground)
    energy_resort = sorted(energy) ##sort energy
    dynamic_energy.append(sorted(energy_resort))
    return energy_resort

def renew_calc_states(filename,filename_new,Nstate):
    regex_1 = re.compile('(?<=root=)[0-9]+')  ##regex assertion (?<=)
    regex = re.compile('td=\(.*root=[0-9]+\)')
    with open (filename,'r') as f:
        with open (filename_new,'w+') as f_new:
            for line in f:
                if not regex.search(line):
                    f_new.write(line)
                else:
                    if Nstate == 0:
                        f_new.write(re.sub(regex,'',line))
                    else:
                        f_new.write(re.sub(regex_1,str(Nstate),line))

def get_element_mass():
    with open('initial_codition', 'r') as f:  # element x y z  p_x p_y p_z
        element_mass = []
        for n, value in enumerate(f):
            data = value.split()[0].capitalize()
            mass_au = float(masses[data] * amu)
            element_mass.append(mass_au)
        return np.array(element_mass)


def get_position_matrix():
    with open('velocity', 'r') as f:
        position_matrix = []
        for n, value in enumerate(f):
            data = list(map(float, value.split()[1:4]))
            data = [i / ang for i in data]
            position_matrix.append(data)
        return np.array(position_matrix)


def get_momentum_matrix():  # element x y z  p_x p_y p_z
    with open('velocity', 'r') as f:
        momentum_matrix = []
        for n, value in enumerate(f):
            data = list(map(float, value.split()[-3:]))
            data = [i / Velo for i in data]
            momentum_matrix.append(data)
        return np.array(momentum_matrix)

"""
Nuclear motion   by velocity verlet  Algorithm
x_n+1 = x_n +  p_n /m_i * Δt  + 1/2mi  * F_n * Δt^2
p_n+1 = p_n +　1/2mi * (F_n+1 +F_n)  * Δt
"""

def update_positon_matrix(position_matrix, momentum_matrix, grad_matrix, element_mass):
    def update_positon(q, p, F, m):
        q_back = q + p / m * step_time + 0.5 * F / m * step_time ** 2
        return q_back
    position_matrix_back = []
    for i in range(natom):
        data = [update_positon(position_matrix[i][j], momentum_matrix[i][j],
                               -grad_matrix[i][j], element_mass[i]) for j in range(3)]
        position_matrix_back.append(data)
    return np.array(position_matrix_back)
## notice atom force symbol postive + or negtive  - 

def update_momentum_matrix(momentum_matrix, grad_matrix, grad_matrix_back, element_mass):
    def update_momentum(p, F, F_b, m):
        p_back = p + 0.50 / m * (F + F_b) * step_time
        return p_back
    velocity_matrix_back = []
    for i in range(natom):
        data = [update_momentum(momentum_matrix[i][j], -grad_matrix[i][j],
                                -grad_matrix_back[i][j], element_mass[i]) for j in range(3)]
        velocity_matrix_back.append(data)
    return np.array(velocity_matrix_back)


def Keyvalue(*args):
    def Dihedral_angle(p1, p2, p3, p4):
        q1 = np.subtract(p2, p1)
        q2 = np.subtract(p3, p2)
        q3 = np.subtract(p4, p3)
        # Calculate cross vectors
        q1_x_q2 = np.cross(q1, q2)
        q2_x_q3 = np.cross(q2, q3)
        # Calculate normal vectors
        n1 = q1_x_q2/np.linalg.norm(q1_x_q2)
        n2 = q2_x_q3/np.linalg.norm(q2_x_q3)
        # Calculate unit vectors
        u1 = n2
        u3 = q2/np.linalg.norm(q2)
        u2 = np.cross(u3, u1)
        # Calculate cosine and sine
        cos_theta = np.dot(n1, u1)
        sin_theta = np.dot(n1, u2)
        theta = -np.arctan2(sin_theta, cos_theta)
        theta_deg = np.degrees(theta)
        return(theta_deg)

    def Bond_angle(p1, p2, p3):
        q1 = np.subtract(p2, p1)
        q2 = np.subtract(p3, p2)
        cos_theta = np.dot(q1, q2) / (np.linalg.norm(q1) * np.linalg.norm(q2))
        theta = np.arccos(cos_theta)
        theta_deg = 180.00 - np.degrees(theta)
        return(theta_deg)

    def Bond_length(p1, p2):
        length = np.linalg.norm(np.subtract(p2, p1))
        return(length)

    N = len([*args])
    n = len(set([tuple(i) for i in [*args]]))  # 变为tuple 可哈希
    if N == n:
        if N == 2:
            return(Bond_length(*args))
        elif N == 3:
            return(Bond_angle(*args))
        elif N == 4:
            return(Dihedral_angle(*args))
        else:
            print('too many paramenters')
            sys.exit()
    else:
        print("atomic coordinates repeat")
        sys.exit()

def on_the_fly():
    current_time =  time.asctime( time.localtime(time.time()) )
    position_matrix_before = [ ]
    momentum_matrix_before = [ ]
    grad_matrix_before = [ ]
    q1 = [ ] ; q2 = [ ] ; q3 = [ ]
    element_mass =  get_element_mass()

    for nloop in range(int(total_time/SI_step_time)):
        if nloop == 0:  #first step 
            position_matrix = get_position_matrix()
            q1 = position_matrix
            software_running()
            print("---------The molecular energy has been calculated at %s  ---------" %(current_time))
            grad_matrix_before =  get_grad_matrix()
            momentum_matrix_before = get_momentum_matrix()
            print( "dynamic program will continue to the next ( %d ) iteration at %s  "%((nloop+1), current_time))
        elif nloop == 1 :  #second step 
            position_matrix = update_positon_matrix(position_matrix_before, momentum_matrix_before, grad_matrix_before, element_mass)
            q1 = q1 
            q2 = position_matrix
            replace_coordinate()
            software_running() 
            print("---------The molecular energy has been calculated completely at %s  ---------" %(current_time))
            grad_matrix =  get_grad_matrix()
            momentum_matrix  = update_momentum_matrix(momentum_matrix_before, grad_matrix, grad_matrix_before, element_mass)
            position_matrix_before = position_matrix  ## This is the p/m/f at the current step (i step ),
            momentum_matrix_before = momentum_matrix  ## used for the next step(i+1 step)  calculation. 
            grad_matrix_before = grad_matrix        ## so the suffix "before" is uesed 
            print( "Dynamic program will continue to the next ( %d ) iteration at %s  " %((nloop+1), current_time))
        else:  #next step 
            position_matrix = update_positon_matrix(position_matrix_before, momentum_matrix_before, grad_matrix_before, element_mass)
            q1 = q2 
            q2 = position_matrix_before
            q3 = position_matrix
            replace_coordinate()
            software_running() 
            print("---------The molecular energy has been calculated completely at %s  ---------" %(current_time))
            force_matrix = get_grad_matrix()
            momentum_matrix = update_momentum_matrix(momentum_matrix_before, force_matrix, grad_matrix_before, element_mass)
            position_matrix_before = position_matrix ## This is the p/m/f at the current step (i step ),
            momentum_matrix_before = momentum_matrix ## used for the next step(i+1 step)  calculation. 
            grad_matrix_before = force_matrix       ## so the suffix "before" is uesed 
            print( "Dynamic program will continue to the next ( %d ) iteration at %s  "%((nloop+1), current_time))
        
def  check_hopping():
    dyn_states = 2 #the current electronic state
    states = np.arange(0,dyn_states)
    hop_type = None
    ##calculate energy difference of the q1 q2 q3 point　 suffix "u" is up  "d" is down 
    if  dyn_states == states[-1] :  ## the highest state  e.g. S1 -> S0
        state_d_num = states -1 
        delta_q1_u = 0.000000 ; delta_q2_u = 0.000000 ; delta_q3_u = 0.000000
        delta_q1_d = (dynamic_energy[nloop -2 ][states] - dynamic_energy[nloop - 2][state_d_num] ) * eV
        delta_q2_d = (dynamic_energy[nloop -1 ][states] - dynamic_energy[nloop - 1][state_d_num] ) * eV
        delta_q3_d = (dynamic_energy[nloop -1 ][states] - dynamic_energy[nloop - 1][state_d_num] ) * eV
    elif dyn_states == states[0]: ##the lowest state  e.g. S0 dynamic 
        state_u_num = 1 
        state_d_num = 0
        delta_q1_d = 0.0000 ;delta_q2_d=0.00000 ; delta_q3_d = 0.00000 
        delta_q1_u = (dynamic_energy[nloop -2 ][1] - dynamic_energy[nloop - 2][0] ) * eV
        delta_q2_u = (dynamic_energy[nloop -1 ][1] - dynamic_energy[nloop - 1][0] ) * eV
        delta_q3_u = (dynamic_energy[nloop -1 ][1] - dynamic_energy[nloop - 1][0] ) * eV
    else :   ##the middle state
        state_u_num = states + 1
        state_d_num = states - 1
        delta_q1_u = (dynamic_energy[nloop -2 ][state_u_num] - dynamic_energy[nloop - 2][state_d_num] ) * eV
        delta_q2_u = (dynamic_energy[nloop -1 ][state_u_num] - dynamic_energy[nloop - 1][state_d_num] ) * eV
        delta_q3_u = (dynamic_energy[nloop -1 ][state_u_num] - dynamic_energy[nloop - 1][state_d_num] ) * eV
        delta_q1_d = (dynamic_energy[nloop -2 ][state_u_num] - dynamic_energy[nloop - 2][state_d_num] ) * eV
        delta_q2_d = (dynamic_energy[nloop -1 ][state_u_num] - dynamic_energy[nloop - 1][state_d_num] ) * eV
        delta_q3_d = (dynamic_energy[nloop -1 ][state_u_num] - dynamic_energy[nloop - 1][state_d_num] ) * eV
    
    """
    check the type of hopping 
    1. A double hop-upward and downward 
    2. A downward hop
    3. A upward hop
    """
    if delta_q1_d > delta_q2_d and delta_q2_d < delta_q3_d :
        if delta_q1_u > delta_q2_u  and delta_q2_u < delta_q3_u and delta_q2_u < kin_energy[nloop - 1 ]:
            if delta_q2_u <= threshold and delta_q2_d <= threshold :
                hop_type = 1 
            elif delta_q2_d <= threshold and delta_q2_u > threshold :
                hop_type = 2 
            elif delta_q2_u <= threshold and delta_q2_d > threshold :
                hop_type = 3 
            else :
                print ("No minmum energy separation smaller than %s ev was found at %d step " %(threshold, nloop))
        else:
            if delta_q2_d <= threshold :
                hop_type = 2 
            else:
                print("The minimum energy separations for downward transitions are larger than  %s" %threshold) 
    elif delta_q1_u > delta_q2_u and delta_q2_u < delta_q3_u and delta_q2_u < kin_energy[nloop - 1]:
        if delta_q2_u <= threshold :
            hop_type = 3 
        else: 
            print (" The minimum energy separations for upward transitions are larger \
                 than %s at %d step" %(threshold ,nloop) ) 
    else: 
        print ("No local minmum between PESs observed at %d step" %nloop)

        
