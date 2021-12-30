import os
import re
import subprocess
import sys
import time

import numpy as np

ang = 0.529177210903
eV = 27.211386245988

__all__ = ["software_running", "read_wavefunction", "delete_wavefunction", "check_initial",
           "get_grad_matrix", "get_energy", "renew_calc_states", "renew_calc_states", "print_traj_cicoe",
           "replace_coordinate"]

try:
    with open("OrcaPath", "r") as f:
        for line in f:
            orca_path = line.split()[0];break
except:
    print("The Orca path of file is not found, try find Orca path from environment(which orca).")
    path = subprocess.run("which orca", shell=True,
                          stdout=subprocess.PIPE, text=True)
    orca_path = path.stdout.replace("\n", "")


def software_running(filename="orca.inp"):
    def current_time():
        return time.asctime(time.localtime(time.time()))
    prefix = filename.split(".")[0]
    output = prefix + ".out"
    print("The ORCA begin running at %s" % current_time(), flush=True)
    file = orca_path + "  " + filename + " &> " + output
    proc = subprocess.Popen(file, shell=True)
    proc.wait()
    print("The ORCA ended at %s" % current_time(), flush=True)


def check_initial(nstates):
    if os.path.isfile("orca.inp"):
        flag_f = False
        flag_r = False
        flag_p = True
        with open("orca.inp", 'r') as f:
            for line in f:
                if re.search("EnGrad", line, re.IGNORECASE):
                    flag_f = True
                if re.search(r'iroot\s+%s' % nstates, line, re.IGNORECASE):
                    flag_r = True
                if re.search("printlevel\s+3", line, re.IGNORECASE):
                    flag_p = True
            if not flag_f:
                raise Exception("Please use key value 'EnGrad'")
            if nstates >= 1:
                if not flag_p:
                    raise Exception("Please use key value 'printlevel 3'")
            if not flag_r:
                raise Exception(
                    "Please check 'iroot %s' in orca.inp" % nstates)
    else:
        print("orca.inp is not found, please check *.inp filename again")
        sys.exit()


def read_wavefunction():
    pass


def delete_wavefunction():
    pass


def get_grad_matrix(filename="orca.out"):
    regex = re.compile("CARTESIAN GRADIENT")
    data = []
    flag_b = False
    flag_g = False
    with open(filename, 'r+') as f:
        for line in f:
            if not flag_g:
                if regex.search(line):
                    flag_g = True
            else:
                if flag_b:
                    if line.strip():
                        data.append(list(map(float, line.split()[3:])))
                    else:
                        break
                else:
                    if not line.strip():
                        flag_b = True
    return np.array(data)


def _get_keyword(filename="orca.inp"):
    regex = re.compile("%TDDFT", re.IGNORECASE)
    key = []
    flag_td = False
    # %TDDFT Nroot 5
    #        iroot 1
    #        tda flase
    #        end
    with open(filename, 'r') as f:
        for line in f:
            if not flag_td:
                if regex.search(line):
                    key.append(line)
                    flag_td = True
            else:
                key.append(line)
                if re.findall("end", line, re.IGNORECASE):
                    break
    return key


word = _get_keyword()


def get_energy(filename='orca.out'):
    regex_scf = re.compile("Total Energy       :")
    regex_l1 = re.compile("lowest eigenvalues of")
    regex_l2 = re.compile("Lowest Energy  ")
    regex_e = re.compile("E\(tot\)  =")
    flag_l1 = False
    flag_l2 = False
    flag_scf = False
    e_ground = 0.0
    all_excited = []
    e_excited = 0.0
    with open(filename, 'r', encoding='UTF-8') as f:
        for line in f:
            if flag_scf:
                if flag_l1:
                    if regex_l2.search(line):
                        flag_l2 = True
                        flag_l1 = False
                        all_excited.append(excited)
                    if not flag_l2:
                        excited.append(float(line.split()[3]) + e_ground)
                else:
                    if regex_l1.search(line):
                        flag_l1 = True
                        flag_l2 = False
                        excited = [e_ground]
            if regex_e.search(line):
                e_excited = float(line.split()[-2])
                break
            else:
                if regex_scf.search(line):
                    e_ground = float(line.split()[3])
                    flag_scf = True
    if all_excited:
        return sorted(all_excited[-1]), e_excited
    else:
        return [e_ground], e_ground


def renew_calc_states(nstates, filename='orca.inp', filename_new=None, st=None, spin=None, charge=None, remove=None, add=None):
    # exclamation(!)!CAM-B3LYP nousesym pal4  D3 RIJCOSX def2-SVP def2/J def2-SVP/C miniprint  TightSCF grid4 EnGrad
    # %pal nprocs 4 end
    # %TDDFT NROOTS 5
    # IROOT    1
    # printlevel 3
    # tda false
    # end
    flag_td = False
    flag_keyword = False
    flag_td_end = False
    flag_file = False
    regex_td = re.compile("%TDDFT", re.IGNORECASE)
    count = 0
    td_n = 0
    with open(filename, 'r') as f:
        for n, line in enumerate(f):
            if regex_td.search(line):
                flag_td = True
                td_n = n
                break
    if not filename_new:
        flag_file = True
        filename_new = "%s.bak" % filename
    with open(filename, 'r') as f, open(filename_new, 'w+') as f_new:
        for n, line in enumerate(f):
            if re.search("^!", line):
                flag_keyword = True
                key_line = n
            else:
                flag_keyword = False
            if flag_keyword:
                if remove:  # remvoe keywords
                    if remove == "force":
                        remove = "EnGrad"
                    if re.search(remove, line, re.IGNORECASE):
                        line = re.sub(remove, "", line)
                if add and count < 1:  # add others keywords
                    if add == "force":
                        add = "EnGrad"
                    line = re.sub('\r?\n', '', line).strip()
                    line = line + " " + add + '\n'
                    count += 1
            else:
                if flag_td:
                    if not flag_td_end:
                        if n >= td_n:
                            if re.search("end", line, re.IGNORECASE):
                                flag_td_end = True
                            if nstates >= 1:
                                # the excited states ->excited states(S2->S1)
                                if re.search("iroot", line, re.IGNORECASE):
                                    line = re.sub(
                                        "iroot\s+[0-9]+", "iroot  %s" % nstates, line, flags=re.IGNORECASE)
                            else:
                                # the excited states -> ground states(S1->S0)
                                line = ""
                else:
                    # ground states -> the excited states
                    if nstates == 0:
                        raise Exception('iroot(%s) must great 0' % nstates)
                    if n == key_line + 1:
                        # write in the after of keyword(!)
                        old_line = "".join(x for x in word)
                        old_line = re.sub(
                            'iroot\s+[0-9]+', "iroot  %s" % nstates, old_line, flags=re.IGNORECASE)
                        line = old_line + line
            f_new.write(line)

    if flag_file:
        os.remove(filename)
        os.rename(filename_new, filename)

# print(get_energy())


def replace_coordinate(new_coordinate, filename='orca.inp', filename_new=None):
    regex = re.compile('[A-Za-z]{1,2}\s*(\s*(-?[0-9]+\.[0-9]*)){3}')
    regex_1 = re.compile('(\s*(-?[0-9]+\.[0-9]*)){3}')
    count = 0
    flag_file = False
    if not filename_new:
        flag_file = True
        filename_new = "%s.bak" % filename
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


def print_traj_cicoe(time, filename='orca.out'):
    with open("traj_cicoe.log", "a+") as f, open(filename, 'r', encoding='UTF-8') as f_o:
        regex = re.compile("the weight of the individual")
        flag = False
        try:
            a = format(time, '>10d')
        except ValueError:
            a = format(time, '>10.2f') + '  fs'
        f.write("simulation time: t=" + a + '\n')
        for line in f_o:
            if not flag:
                if regex.search(line):
                    flag = True
            else:
                if '-------' in line:
                    break
                else:
                    f.write(line)
