import os
import re
import sys
import time

import numpy as np

eV = 27.21138602
ang = 0.529177257507
numlti = {'Singlet': 1,
          'Doublet': 2,
          'Triplet': 3,
          'Quartet': 4,
          'Quintet': 5,
          'Sextet': 6,
          'Septet': 7,
          'Octet': 8
          }


def current_time():
    return time.asctime(time.localtime(time.time()))


def software_running(filename='gauss.gjf'):
    print("The Gaussian begin running at %s" % current_time(), flush=True)
    file = 'g16  ' + filename
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
                        line + ' '+'geom=allcheck ' + '\n')
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
                if re.search("root", line, re.IGNORECASE):
                    flag_r = True
            if not flag_p:
                print("Plase use detailed output(#P) ")
                sys.exit()
            if not flag_f:
                print("Please use key value 'force'")
                sys.exit()
        #if  states_involved >= 2:
        #    if not flag_r:
        #        print("Please use key value 'root'")
        #        sys.exit()
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


def get_grad_matrix(filename, natom):
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
                        data.append(list(map(float, line.split()[2:])))
                        if len(data) == natom:  # mathch line
                            break
            else:
                addflag = True
    # the gradient is opposite to the direction of force
    gradient_matrix = -np.array(data)
    return gradient_matrix


def get_energy(filename):
    """
    get S0/the excited stated energy from detailed output file. 
    """
    regex = re.compile('SCF Done')
    regex_1 = re.compile('Excitation Energies \[eV] at current iteration:')
    regex_2 = re.compile('Convergence achieved on expansion vectors')
    regex_3 = re.compile('Total Energy')
    #regex_4 = re.compile('?<=(Excited State(\s{3}[1-9]|\s{2}1[0-9]):\s{6})(Triplet|Singlet)')
    #regex_4 = re.compile('Excited State\s*[0-9]*\s*:')
    energy = []
    tmp = []
    E_ground = 0.0
    E_total = 0.0
    flag_g = True
    flag_e = False
    flag_c = False
    #nmulti = [ ]
    #nstate = [ ]
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
                    # if regex_4.search(line):
                     #   nmulti.append(line.split()[3][:-2])
                else:
                    if flag_e:
                        if regex_2.search(line):
                            flag_c = True
                            for i in tmp:
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


def renew_calc_states(nstate, filename, filename_new=None):
    """
    This is the function to adjust the interested state
    e.g. S2-S1 ,S1-S0, S0-S2
    when the input file is the ground state(S0),
    it will write the previous keywords to calculate the excited state energy
    """
    # regex assertion (?<=)
    regex_1 = re.compile('(?<=root=)[0-9]+', re.IGNORECASE)
    regex = re.compile('(tda?=?\(.*?\)|tda)', re.IGNORECASE)
    flag = False
    flag_td = False
    with open(filename) as f:
        for line in f:
            if regex.search(line):
                flag_td = True
    if not filename_new:
        flag = True
        filename_new = "%s.bak" % filename
    if flag_td:
        with open(filename, 'r') as f, open(filename_new, 'w+') as f_new:
            for line in f:
                if not regex.search(line):
                    f_new.write(line)
                else:
                    if nstate == 0:
                        f_new.write(re.sub(regex, '', line))
                    else:
                        f_new.write(re.sub(regex_1, str(nstate), line))
    else:  # S0-> the excited states
        regex_2 = re.compile('force', re.IGNORECASE)
        with open(filename, 'r') as f,  open(filename_new, 'w+') as f_new:
            for line in f:
                if not regex_2.search(line):
                    f_new.write(line)
                else:
                    line = re.sub('\r?\n', '', line).strip()
                    WORD = re.sub(regex_1, '1', word)
                    f_new.write(line + ' '+WORD + '\n')
    if flag:
        os.remove(filename)
        os.rename(filename_new, filename)


def replace_coordinate(new_coordinate, filename='gauss.gjf', ):
    regex = re.compile('[A-Za-z]{1,2}\s*(\s*(-?[0-9]+\.[0-9]*)){3}')
    regex_1 = re.compile('(\s*(-?[0-9]+\.[0-9]*)){3}')
    count = 0
    with open(filename, 'r') as f, open("%s.bak" % filename, 'w+') as f_new:
        for n, line in enumerate(f):
            if not regex.findall(line):
                f_new.write(line)
            else:
                replace_coord = ''.join(format(i * ang, '>18.10f')
                                        for i in new_coordinate[count][:3])
                f_new.write(re.sub(regex_1, replace_coord, line))
                count += 1
    os.remove(filename)
    os.rename("%s.bak" % filename, filename)


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


def print_traj_cicoe(time,filename='gauss.log'):
    with open('traj_cicoe.log', 'a+') as f:
        regex = re.compile('Excited State\s*[0-9]*:')
        regex_1 = re.compile('SavETr')
        f.write('simulation time: t=' +
                format(time, '>10.2f') + '\n')
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
