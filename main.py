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
import os
import shutil
import time

import numpy as np

from interface import *
# from interface import eV, ang
from interface.prepare import SI_step_time, total_time, dynamics_time, nloop
from interface.prepare import dyn_states, states_min, states_max, states_involved
from interface.prepare import element, element_mass, atom_list
from interface.prepare import inputfile, inputfile_d, inputfile_u, outputfile_d, outputfile_u
from interface.prepare import natom, initial_mom, initial_coord
from interface.prepare import soft, threshold
from interface.prepare import wfnfile, wfn_suffix, EnterDir
from ultis.geom import coord_corrections, coord_mom_rot, coord_corrections, mom_corrections
from ultis.geom import critical_value, calculate_kinetic
from ultis.velocity_verlet import update_momentum_matrix, update_position_matrix

# All units in the code are a.u./hartree
# The units of initial coordination/velocity is Angstrom/Å  and Bohr/a.u.
np.set_printoptions(precision=8)
potential_energy = []  # molecular dynamics energy
kin_energy = []  # molecular  kinetic energy
total_energy = []  # kinetic + potential
q1 = []
q2 = []
q3 = []
coord = []  # x y z
momentum = []  # px py pz
rot_matrix = []
grad = []  #
scf_energy = []
hop_time = [0]  # hop time in dynamic programs


# flag_gap = False


# software 1.GAUSSIAN  2.ORCA  3.BDF 4.MOLPRO 5.MOLCAS
# soft = software


def current_time():
    return time.asctime(time.localtime(time.time()))


def get_initial_condition():
    coord.append(initial_coord)
    momentum.append(initial_mom)


def analyse_result():
    """
    grad, potential_energy, total_energy
    """
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
        renew_calc_states(states_involved - 1, inputfile,
                          remove="force")  # root=states_involved - 1
        software_running()
        E_exc, E_scf = get_energy()
        # get the excited stated excited energy
        potential_energy.append(E_exc)
        delete_wavefunction()
        #  renew_calc_states(dyn_states, inputfile) # error  *****
        renew_calc_states(dyn_states, inputfile, add="force")


def print_matrix(value):
    for i in range(natom):
        ele = format(element[i], '<5s')
        xyz = "".join(format(x, '>18.10f') for x in value[i])
        print(ele + xyz, flush=True)


def on_the_fly():  # 在software check——hooping 之后
    """
    renew the current step momentum and the nex step position
    Nuclear motion by velocity verlet Algorithm
    x_n+1 = x_n + p_n /m_i * Δt  + 1/2mi * F_n * Δt^2
    p_n+1 = p_n + 0.5 * (F_n+1 +F_n)  * Δt
    """
    global q1, q2, q3
    if nloop == 0:  # first step
        # coord.append(update_position_matrix(coord[nloop], momentum[nloop], grad[nloop], element_mass))
        # 2step position
        next_coord = update_position_matrix(coord[nloop], momentum[nloop], grad[nloop])
        next_coord_c, rot = coord_corrections(coord[0], next_coord, element_mass)
        coord.append(next_coord_c)
        rot_matrix.append(rot)
        q1 = coord[nloop]
        q2 = coord[nloop + 1]
    else:  # next step
        # the next position  & the current step momentum
        # momentum.append(update_momentum_matrix(momentum[nloop - 1], grad[nloop], grad[nloop - 1]))
        # coord.append(update_position_matrix(coord[nloop], momentum[nloop], grad[nloop]))
        next_mom = update_momentum_matrix(momentum[nloop - 1], grad[nloop], grad[nloop - 1])
        next_mom_c = next_mom @ rot_matrix[nloop - 1]
        momentum.append(next_mom_c)
        next_coord = update_position_matrix(coord[nloop], momentum[nloop], grad[nloop])
        next_coord_c, rot = coord_corrections(coord[0], next_coord, element_mass)
        # next_mom_c = mom_corrections(momentum[0], next_mom)
        # print("This difference between next_mom and next_mom_c is %14.8f - %14.8f = %14.8f " % (
        # calculate_kinetic(next_mom), calculate_kinetic(next_mom_c),
        # calculate_kinetic(next_mom) - calculate_kinetic(next_mom_c)))
        coord.append(next_coord_c)
        rot_matrix.append(rot)
        q1 = coord[nloop - 1]
        q2 = coord[nloop]
        q3 = coord[nloop + 1]


def save_wavefunction():
    with EnterDir(soft + "Temp"):
        c = wfn_suffix[soft]
        if nloop == 0:
            shutil.copy(wfnfile, "q1." + c)
        elif nloop == 1:
            shutil.copy(wfnfile, "q2." + c)
        elif nloop == 2:
            shutil.copy(wfnfile, "q3." + c)
        else:
            shutil.copy("q2." + c, "q1." + c)
            shutil.copy("q3." + c, "q2." + c)
            shutil.copy(wfnfile, "q3." + c)


def check_hopping():  # when nloop >=2,begin check hopping
    # global dyn_states
    # global nloop
    hop_type = None
    states = np.arange(states_involved)  # (0,1)
    # calculate energy difference of the q1 q2 q3 point　 suffix "u" is up "d" is down

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
        print("The current states(%d) is the lowest state" %
              dyn_states, flush=True)
        if not (dyn_states == 0 and states_involved <= 1):  # S0 dynamic after the hopping
            print(
                "The energy differences between %s and %s states at the least three point(q1,q2,q3-current) are: "
                "%12.6f %12.6f %12.6f "
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
        with EnterDir(soft + "Temp"):
            c = wfn_suffix[soft]
            # read q1 wavefunction file
            shutil.copy("q1." + c, wfnfile)
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
        with EnterDir(soft + "Temp"):
            c = wfn_suffix[soft]
            # read q1 wavefunction file
            shutil.copy("q1." + c, wfnfile)
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

        # calculate q3, q1 and q2 up state grad

        hop_direction = 'U'
        # hop_p = 0

        # calculate = q3 up_state grad
        renew_calc_states((dyn_states + 1), inputfile, inputfile_u)
        software_running(inputfile_u)
        grad_q3u = get_grad_matrix(outputfile_u)

        # calculate = q1 up_state grad
        with EnterDir(soft + "Temp"):
            c = wfn_suffix[soft]
            # read q1 wavefunction file
            shutil.copy("q1." + c, wfnfile)
        replace_coordinate(q1)
        renew_calc_states((dyn_states + 1), inputfile, inputfile_u)
        software_running(inputfile_u)
        grad_q1u = get_grad_matrix(outputfile_u)

        hop_p_u, delta_grad_q2_u, mom_direction_factor_u, grad_q2m, grad_q2u = get_hop_factor(
            grad_q3m, grad_q3u, grad_q1m, grad_q1u, hop_direction)

        # calculate q3, q1 and q2 down state grad

        hop_direction = 'D'
        # calculate = q3 down_state grad
        with EnterDir(soft + "Temp"):
            c = wfn_suffix[soft]
            # read q1 wavefunction file
            shutil.copy("q3." + c, wfnfile)
        replace_coordinate(q3)
        renew_calc_states((dyn_states - 1), inputfile, inputfile_d)
        software_running(inputfile_d)
        grad_q3d = get_grad_matrix(outputfile_d)

        # calculate = q1 down_state grad
        with EnterDir(soft + "Temp"):
            c = wfn_suffix[soft]
            # read q1 wavefunction file
            shutil.copy("q1." + c, wfnfile)
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
    Ex = 0.0
    Vx = 0.0
    if hop_direction == 'D':  # downward  up
        Ex = (potential_energy[nloop][dyn_states] + potential_energy[nloop][dyn_states - 1]) * 0.50000
        Vx = (potential_energy[nloop][dyn_states] - potential_energy[nloop][dyn_states - 1]) * 0.50000
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
        Ex = (potential_energy[nloop][dyn_states] + potential_energy[nloop][dyn_states + 1]) * 0.50000
        Vx = (potential_energy[nloop][dyn_states + 1] - potential_energy[nloop][dyn_states]) * 0.50000
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
            delta_grad_q2 += (grad_q2u[i][j] - grad_q2d[i][j]) ** 2 / element_mass[i]
            # \frac{F^2_i - F^1_1 }{\sqrt{m_i}}
            mom_direction_factor[i][j] = (grad_q2u[i][j] - grad_q2d[i][j]) / np.sqrt(element_mass[i])
            # \sum{F^1_i * F^2_i}
            F12 += grad_q2d[i][j] * grad_q2u[i][j]

    # calculate parallel/vertical kinetic
    hop_direction_normal = np.zeros((natom, 3))
    hop_direction_matrix = mom_direction_factor / np.sqrt(delta_grad_q2)  # s_i
    mom_parallel = np.zeros((natom, 3))
    kinetic_parallel = 0.0
    # before the hopping
    for i in range(natom):
        # n_i Normalized Si
        hop_direction_normal[i] = hop_direction_matrix[i] / np.linalg.norm(hop_direction_matrix[i])
        # p_parallel = (n_i . p_i) * n_i
        mom_parallel[i] = np.dot(hop_direction_normal[i], momentum[nloop - 1][i]) * hop_direction_normal[i]
        # Parallel kinetic energy
        kinetic_parallel += 0.50 * np.linalg.norm(mom_parallel[i]) ** 2 / element_mass[i]
        # kinetic_total += 0.50 * np.linalg.norm(momentum[nloop - 1]) ** 2 / element_mass[i]

    kinetic_total = calculate_kinetic(momentum[nloop - 1])
    kinetic_vertical = kinetic_total - kinetic_parallel

    print("The total kinetic, parallel kinetic, vertical kinetic is %12.8f, %12.8f , %12.8f (eV)" % (
        kinetic_total * eV, kinetic_parallel * eV, kinetic_vertical * eV), flush=True)

    # calculate a^2 and b^2
    f_aa = delta_grad_q2 / (16 * Vx ** 3)  # a.u. Reduced Planck constant = 1
    f_bb = (scf_energy[nloop - 1] + kinetic_parallel - Ex) / (2 * Vx)
    # notice "E =  quantum energy + parallel kinetic"

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
    if hop_direction == "D":
        print("The current(upper) states %s " % dyn_states, flush=True)
        print("The downer states %s" % (dyn_states - 1), flush=True)
    if hop_direction == "U":
        print("The current(downer) states %s" % dyn_states, flush=True)
        print("The upper states is states %s" % (dyn_states + 1), flush=True)
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
    print("The (F2-F1)^2/m ,F12 is %18.10f %18.10f" % (delta_grad_q2, F12), flush=True)
    print("The energy E_p , Vx, Ex is %18.10f %18.10f %18.10f" % (scf_energy[nloop - 1] + kinetic_parallel, Vx, Ex),
          flush=True)
    print("a^2, b^2, hop_p is %18.10f %18.10f %18.10f" % (f_aa, f_bb, hop_p), flush=True)

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
        raise ValueError("The hop direction is error")

    # delta_energy_q2 = potential_energy_q2(+) - potential_energy_q2(-)
    delta_energy_q2 = potential_energy[nloop - 1][dyn_states] - potential_energy[nloop - 1][last_states]

    # before the hopping
    mom_parallel = np.zeros((natom, 3))
    mom_renew = np.zeros((natom, 3))
    for i in range(natom):
        # n_i Normalized Si
        hop_direction_normal[i] = hop_direction_matrix[i] / np.linalg.norm(hop_direction_matrix[i])
        # p_parallel = (n_i . p_i) * n_i
        mom_parallel[i] = np.dot(hop_direction_normal[i], momentum[nloop - 1][i]) * hop_direction_normal[i]
        # Parallel kinetic energy
        kinetic_parallel += 0.50 * np.linalg.norm(mom_parallel[i]) ** 2 / element_mass[i]

    if hop_direction == 'D':
        # higher state : e.g. S1->S0 k > 1
        increment_k = np.sqrt(1 + delta_energy_q2 / kinetic_parallel)
    else:
        # hop_direction == 'U': Lower state - higher state e.g. :S1->S2 k < 1
        increment_k = np.sqrt(1 - delta_energy_q2 / kinetic_parallel)
    # calculate the momentum after the hopping
    for i in range(natom):
        mom_renew[i] = momentum[nloop - 1][i] + (increment_k - 1) * \
                       np.dot(momentum[nloop - 1][i], hop_direction_normal[i]) * hop_direction_normal[i]

    print("The increment factor(k) for momentum adjustment is %.8f" %
          increment_k, flush=True)
    print("***************")
    # print("Normalized momentum direction factor(S_i) is ", flush=True)
    # print_matrix(mom_direction_factor)
    print("Normalized momentum direction factor(S_i) ", flush=True)
    print_matrix(hop_direction_matrix)
    print("Denormalized momentum direction factor (n_i)", flush=True)
    print_matrix(hop_direction_normal)
    print("The q2 old momentum is", flush=True)
    print_matrix(momentum[nloop - 1])
    print("The parallel momentum is %.8f eV" %
          (kinetic_parallel * eV), flush=True)
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

    # check q2 kinetic energy before and after hopping
    q2_old_k = calculate_kinetic(momentum[nloop - 1])
    # delete q2 momentum and renew q2 momentum
    momentum.pop()
    mom_renew_c = mom_renew @ rot_matrix[nloop - 1]  # ??????
    momentum.append(mom_renew_c)
    q2_new_k = calculate_kinetic(momentum[nloop - 1])

    print("The parallel kinetic is %.8f eV" %
          (kinetic_parallel * eV), flush=True)
    print("The q2 old/new kinetic energy is %15.10f eV %15.10f eV" %
          (q2_old_k * eV, q2_new_k * eV), flush=True)
    print("The kinetic energy difference is %15.8f eV" %
          ((q2_old_k - q2_new_k) * eV), flush=True)
    print("The energy difference is %15.10f eV" %
          (delta_energy_q2 * eV), flush=True)
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

    # delete q3 scf_energy
    scf_energy.pop()  # ??????

    # delete q3 rot_matrix
    rot_matrix.pop()

    # revise q2 kinetic-energy potential-energy total energy
    # delete q2 coordinate and momentum
    with open('simulation.xyz', 'r') as f, open("%s.bak" % 'simulation.xyz', 'w+') as f_new:
        lines = f.readlines()
        f_new.writelines(lines[:-(natom + 2)])
        f_new.write(str(natom) + '\n')
        f_new.write('simulation time: t  =' + format(dynamics_time - SI_step_time, '>10.2f') + ' \n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            xyz = ''.join(format(x * ang, '>18.10f') for x in coord[nloop - 1][i])
            mom = ''.join(format(x, '>18.10f') for x in momentum[nloop - 1][i])
            f_new.write(ele + xyz + mom + '\n')
    os.remove('simulation.xyz')
    os.rename("%s.bak" % 'simulation.xyz', 'simulation.xyz')

    # delete q2 old gradient and write q2 new old gradient in "traj_grad.log"
    with open('traj_grad.log', 'r') as f, open("%s.bak" % 'traj_grad.log', 'w+') as f_new:
        lines = f.readlines()
        f_new.writelines(lines[:-(natom + 2)])
        f_new.write(str(natom) + '\n')
        f_new.write('simulation time: t  =' +
                    format(dynamics_time - SI_step_time, '>10.2f') + ' \n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            gradient = ''.join(format(x, '>18.10f') for x in new_grad_q2[i])
            f_new.write(ele + gradient + '\n')
    os.remove('traj_grad.log')
    os.rename("%s.bak" % 'traj_grad.log', 'traj_grad.log')

    # record hop time/coord/momentum (q2, nloop -1)
    with open('traj_hopping.log', 'a+') as f:
        f.write('{:<8.2f}fs{:>5s} => {:<5s}'.format(dynamics_time - SI_step_time,
                                                    str(before_dyn_states), str(dyn_states)) + '\n')
        f.write('{:<15.8f}{:>15.8f}{:>15.8f}{:>15.8f}'.format(total_energy[nloop - 1],
                                                              potential_energy[nloop - 1][before_dyn_states],
                                                              potential_energy[nloop - 1][dyn_states],
                                                              kin_energy[nloop - 1]))
        f.write('   ' + '\n')
    hop_xyz = 'hop-coord-' + format(dynamics_time - SI_step_time, '.1f') + '.xyz'
    with open(hop_xyz, 'w+') as f:
        f.write(str(natom) + "\n")
        f.write('Hopping Time = ' +
                format(dynamics_time - SI_step_time, '>10.2f') + '\n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            xyz = ''.join(format(x * ang, '>18.10f') for x in coord[nloop - 1][i])
            mom = ''.join(format(x, '>18.10f') for x in mom_renew_c[i])
            f.write(ele + xyz + mom + '\n')

    # revise q2 traj_energy.log (q2) q3 is current coordinate
    with open('traj_energy.log', 'r') as f, open("%s.bak" % 'traj_energy.log', 'w+') as f_new:
        lines = f.readlines()
        f_new.writelines(lines[:-1])
        energy_involved = ''.join(
            format(potential_energy[nloop - 1][i], '>18.8f') for i in range(states_involved))
        test = '{:<8.2f}{:>18.8f}{:>18.8f}{:>18.8f}'.format(dynamics_time - SI_step_time, total_energy[nloop - 1],
                                                            potential_energy[nloop - 1][dyn_states],
                                                            kin_energy[nloop - 1])
        f_new.write(test + '     ' + energy_involved + '\n')
    os.remove('traj_energy.log')
    os.rename("%s.bak" % 'traj_energy.log', 'traj_energy.log')

    # recalculate q3 gradient potential_energy gradient q2 coordinate, gradient momentum
    # coord.append(update_position_matrix(coord[nloop - 1], momentum[nloop - 1], grad[nloop - 1], element_mass))
    # coord.append(update_position_matrix(q2, mom_renew, new_grad_q2, element_mass))
    global q3
    q3 = update_position_matrix(q2, mom_renew, new_grad_q2)
    q3_c, rot = coord_corrections(coord[0], q3, element_mass)
    coord.append(q3_c)
    rot_matrix.append(rot)
    replace_coordinate(coord[nloop])
    software_running()
    save_wavefunction()
    analyse_result()

    # q3_kin = calculate_kinetic(update_momentum_matrix(momentum[nloop - 1], grad[nloop], new_grad_q2))
    q3_kin = calculate_kinetic(update_momentum_matrix(mom_renew, grad[nloop], new_grad_q2))
    q3_scf = scf_energy[nloop]
    q3_total = q3_kin + q3_scf

    print("The q3 kinetic/scf energy is %18.8f eV %18.10f" % (q3_kin * eV, q3_scf), flush=True)
    print("The q3 total energy is %18.10f" % q3_total, flush=True)
    print("*************", flush=True)
    print("The q3 old coordinate(a.u.) is ", flush=True)
    print_matrix(q3_old)
    print("The q3 adjusting coordinate(a.u.) is", flush=True)
    print_matrix(coord[nloop])
    print("The q3 grad is ", flush=True)
    print_matrix(grad[nloop])
    print("*************", flush=True)


# def check_result():
#     # calculate bond length angle  dihedral
#     # global SI_step_time, step_time, flag_gap
#     bond1 = [25, 28]
#     bond2 = [22, 28]
#     bond3 = [22, 31]
#     bond4 = [31, 30]
#     key1 = [coord[nloop][i - 1] * ang for i in bond1]
#     key2 = [coord[nloop][i - 1] * ang for i in bond2]
#     key3 = [coord[nloop][i - 1] * ang for i in bond3]
#     key4 = [coord[nloop][i - 1] * ang for i in bond4]
#     OH1 = critical_value(*key1)
#     OH2 = critical_value(*key2)
#     OH3 = critical_value(*key3)
#     OH4 = critical_value(*key4)
#     if not flag_gap:
#         if np.abs(OH1 - OH2) <= 0.100 or np.abs(OH3 - OH4) <= 0.100:
#             SI_step_time = 0.2
#             step_time = SI_step_time / fs
#             flag_gap = True
#     else:
#         if np.abs(OH1 - OH2) > 0.100 and np.abs(OH3 - OH4) > 0.100:
#             SI_step_time = initial_step_time
#             step_time = SI_step_time / fs
#             flag_gap = False


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
        f.write('simulation time: t  =' + format(dynamics_time, '>10.2f') + ' \n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            xyz = ''.join(format(x * ang, '>18.10f') for x in coord[nloop][i])
            mom = ''.join(format(x, '>18.10f') for x in momentum[nloop][i])
            f.write(ele + xyz + mom + '\n')

    # print traj_grad.log
    with open("traj_grad.log", 'a+') as f:
        f.write(str(natom) + '\n')
        f.write('simulation time: t  =' +
                format(dynamics_time, '>10.2f') + ' \n')
        for i in range(natom):
            ele = format(element[i], '<5s')
            gradient = ''.join(format(x, '>18.10f') for x in grad[nloop][i])
            f.write(ele + gradient + '\n')

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
    # basis_constant()
    # set_var()
    # check_initial(dyn_states)
    get_initial_condition()  # read coordinates from initial condition
    replace_coordinate(coord[0])
    global dynamics_time, nloop
    while dynamics_time <= total_time:
        print("Trajectory calculation of step %s loop start at the %5.2f fs (%s)" % (
            nloop, dynamics_time, current_time()), flush=True)
        software_running()
        save_wavefunction()
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