# This is some function which could deal with GAUSSIAN16 
# eg: run GAUSSIAN16 software： software_running()
# eg: get gradient matrix (3 * N_atom): get_grad_matrix()
# eg: change keyword in the GAUSSIAN input file : renew_calc_states()
# eg: replace_coordinate()
# eg: print_traj_cicoe()
# Author: zibo wu <zbwu1996@gmail.com>

import re
import sys
import time
import subprocess
import numpy as np
import os 

ang = 0.529177210903
eV = 27.211386245988

__all__ = ["software_running","read_wavefunction", "delete_wavefunction", "check_initial",
            "get_grad_matrix", "replace_coordinate", "renew_calc_states", "analyse_result"
            "print_traj_cicoe"]

def software_running(filename='gauss.gjf'):
    def current_time():
        return time.asctime(time.localtime(time.time()))
    print("The Gaussian begin running at %s" % current_time(), flush=True)
    file = 'g16 ' + filename
    proc = subprocess.Popen(file, shell=True)
    proc.wait()
    print("The Gaussian ended at %s" % current_time(), flush=True)


def read_wavefunction(filename='gauss.gjf'):
    flag_chk = False
    with open(filename, 'r') as f:
        for line in f:
            if re.search('geom=allcheck', line, re.IGNORECASE):
                flag_chk = True
    regex = re.compile("force", re.IGNORECASE)
    with open(filename, 'r') as f, open("%s.bak" % filename, 'w+') as f_new:
        for line in f:
            if not regex.search(line):
                f_new.write(line)
            else:
                if not flag_chk:
                    line = re.sub('\r?\n', '', line).strip()
                    f_new.write(
                        line + ' ' + 'geom=allcheck ' + '\n')
                else:
                    f_new.write(line)
    os.remove(filename)
    os.rename("%s.bak" % filename, filename)


def delete_wavefunction(filename='gauss.gjf'):
    regex = re.compile('geom=allcheck', re.IGNORECASE)
    with open(filename, 'r') as f, open("%s.bak" % filename, 'w+') as f_new:
        for line in f:
            if not regex.search(line):
                f_new.write(line)
            else:
                f_new.write(re.sub(regex, '', line))
    os.remove(filename)
    os.rename("%s.bak" % filename, filename)


def check_initial():
    if os.path.isfile("gauss.gjf"):
        flag_f = False
        flag_p = False
        flag_r = False
        with open("gauss.gjf", 'r') as f:
            for line in f:
                if re.match("#([pP])", line):
                    flag_p = True
                if re.search("force", line, re.IGNORECASE):
                    flag_f = True
                # e.search("root", line, re.IGNORECASE):
                # flag_r = True
            if not flag_p:
                print("Please use detailed output(#P) ")
                sys.exit()
            if not flag_f:
                print("Please use key value 'force'")
                sys.exit()
    else:
        print("gauss.gjf is not found, please check *.gjf filename again")
        sys.exit()
    if os.path.isfile("initial_condition"):
        with open("initial_condition", 'r') as f:
            for line in f:
                if len(line.split()) == 0:
                    continue
                else:
                    if len(line.split()) != 7:
                        print(
                            "Please check the format of 'initial_condition' (elem x y z px py pz)")
                        sys.exit()
    else:
        print("'initial_condition' is not found ,please check it again")


def _get_keyword(filename='gauss.gjf'):
    regex = re.compile('(tda?=?\(.*?\)|tda)', re.IGNORECASE)
    key = []
    with open(filename) as f:
        for line in f:
            if regex.search(line):
                keyword = ' '.join(str(x) for x in regex.findall(line))
                key.append(keyword)
    tda = ' '.join(str(x) for x in key)
    return tda.strip()


word = _get_keyword()


def get_grad_matrix(filename='gauss.log'):
    regex = re.compile('Forces \(Hartrees/Bohr\)')
    data = []
    addflag = False
    count = 0
    with open(filename, 'r') as f:
        for line in f:
            if not regex.search(line) or addflag:
                if addflag:
                    if count < 2:
                        count += 1
                    else:
                        if '---' in line:
                            break
                        else:
                            data.append(list(map(float, line.split()[2:])))
            else:
                addflag = True
    # the gradient is opposite to the direction of force
    gradient_matrix = -np.array(data)
    return gradient_matrix


def get_energy(filename='gauss.log'):
    """
    get S0/the excited stated energy from detailed output file. 
    """
    regex = re.compile('SCF Done')
    regex_1 = re.compile('Excitation Energies \[eV] at current iteration:')
    regex_2 = re.compile('Convergence achieved on expansion vectors')
    regex_3 = re.compile('Total Energy')
    energy = []
    tmp = []
    E_ground = 0.0
    E_total = 0.0
    flag_g = True
    flag_e = False
    flag_c = False
    with open(filename, 'r') as f:
        for line in f:
            if flag_g:
                if regex.search(line):  # SCF Done
                    E_ground = float(line.split()[4])
                    energy.append(E_ground)
                    flag_g = False
            else:
                if flag_c:  # excited states iteration done
                    if regex_3.search(line):
                        E_total = float(line.split()[4])
                else:
                    if flag_e:
                        if regex_2.search(line):
                            flag_c = True
                            for i in tmp:
                                # get every states energy/eV
                                energy.append(
                                    float(i.split()[3]) / eV + E_ground)
                        else:
                            if regex_1.search(line) and tmp:
                                tmp = []
                            else:
                                tmp.append(line)
                    else:
                        if regex_1.search(line):
                            flag_e = True
        if len(energy) == 1:
            E_total = E_ground
    return sorted(energy), E_total


def renew_calc_states(nstate, filename='gauss.gjf', filename_new=None, st=None, spin=None, charge=None, remove=None, add=None):
    """
    This is the function to adjust the interested state,spin-multiplicity and charge
    e.g. S2-S1 ,S1-S0, S0-S2, S0-T1
    when the input file is the ground state(S0),
    it will read the previous keywords to calculate the excited state energy
    """
    regex_td = re.compile('(tda?=?\(.*?\)|tda)', re.IGNORECASE)
    regex_root = re.compile('(?<=root=)[0-9]+', re.IGNORECASE)
    regex_st = re.compile("singlets|triplets|50-50", re.IGNORECASE)
    regex_spin = re.compile("^\s*?-?[0-9]\s+[1-9]\s*$")
    flag_td = False
    flag_keyword = False
    flag_file = False
    ST = {'S': 'singlets', 'T': 'triplets', 'ST': '50-50'}
    with open(filename) as f:
        for line in f:
            if regex_td.search(line):
                flag_td = True
                break
    if not filename_new:
        flag_file = True
        filename_new = "%s.bak" % filename
    with open(filename, 'r') as f, open(filename_new, 'w+') as f_new:
        for line in f:
            if re.search('#(p|P|s|S)?', line):
                flag_keyword = True
            if line.isspace():
                if flag_keyword:
                    flag_keyword = False
            if flag_keyword: #change the key word
                if flag_td:
                    if nstate == 0:  # excited states -> ground states
                        line = re.sub(regex_td, '', line)
                    else:
                        line = re.sub(regex_root, str(nstate), line)
                else:  # ground states -> the excited states
                    if re.search('#(p|P|s|S)?', line):
                        if nstate >= 1:
                            line = re.sub('\r?\n', '', line).strip()
                            WORD = re.sub(regex_root, str(nstate), word)
                            line = line + ' ' + WORD + '\n'
                if st:  # change the singlet/Triplets excited states for closed-shell systems
                    if re.search(regex_td, line):
                        if re.search(regex_st, line):
                            line = re.sub(regex_st, ST[st], line)
                        else:
                            LINE = re.search(
                                "tda?=\(", line).group() + ST[st] + ","
                            line = re.sub("tda?=\(", LINE, line, re.IGNORECASE)
                if remove:
                    if re.search(remove, line, re.IGNORECASE):
                        line = re.sub(remove, line)
                if add:
                    if re.search('#(p|P|s|S)?', line):
                        line = re.sub('\r?\n', '', line).strip()
                        line = line + " " + add + '\n'
            else:
                if regex_spin.search(line):
                    c, s = [int(i) for i in line.split()]
                    if not charge:
                        charge = c
                    if not spin:
                        spin = s
                    line = str(charge) + " " + str(spin) + "\n"
            f_new.write(line)
    if flag_file:
        os.remove(filename)
        os.rename(filename_new, filename)


def replace_coordinate(new_coordinate, filename='gauss.gjf',filename_new=None):
    regex = re.compile('[A-Za-z]{1,2}\s*(\s*(-?[0-9]+\.[0-9]*)){3}')
    regex_1 = re.compile('(\s*(-?[0-9]+\.[0-9]*)){3}')
    count = 0
    flag_file = False
    if not filename_new:
        flag_file = True
        filename_new = "%s.bak" %filename
    with open(filename, 'r') as f, open(filename_new, 'w+') as f_new:
        for n, line in enumerate(f):
            if not regex.findall(line):
                f_new.write(line)
            else:
                replace_coord = ''.join(format(i * ang, '>18.10f')
                                        for i in new_coordinate[count][:3])
                f_new.write(re.sub(regex_1, replace_coord, line))
                count += 1
    if flag_file:
        os.remove(filename)
        os.rename(filename_new, filename)
        


def analyse_result(filename='gauss.log'):
    regex = re.compile('Error termination')
    regex_1 = re.compile('Excited State\s*[0-9]*:')
    flag = False
    with open(filename, 'r') as f:
        for line in f:
            if regex.search(line):
                print('Gaussian possible convergence error', flag=True)
                flag = True
            if regex_1.search(line):
                if float(line.split()[4]) < 0:
                    print('Excitation energy is negative', flush=True)
                    flag = True
    return flag


def print_traj_cicoe(time, filename='gauss.log', flag: bool = True):
    with open('traj_cicoe.log', 'a+') as f:
        regex = re.compile('Excited State\s*[0-9]*:')
        regex_1 = re.compile('SavETr')
        if flag:
            f.write('simulation time: t=' + format(time, '>10.2f') + '\n')
        else:
            f.write('simulation time: t=' + format(time, '>10d') + '\n')
        f.write('\n')
        addflag = False
        with open(filename, 'r+') as f_1:
            for line in f_1:
                if not regex.search(line) or addflag:
                    if addflag:
                        if regex_1.search(line):
                            f.write('\n')
                            break
                        else:
                            f.write(line)
                else:
                    addflag = True
                    f.write(line)
