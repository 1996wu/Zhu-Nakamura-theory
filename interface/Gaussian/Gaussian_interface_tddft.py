import os
import re
import numpy as np
from abc import ABC, abstractmethod

from ..prepare import soft
from ..prepare import output_suffix, input_suffix, input_prefix
from ..PublicFunction import PublicFunction

eV = 27.21138602
ang = 0.529177257507  # Angstrom/Å e-10


class Gaussian(PublicFunction):
    def energy(self, filename, nstates: int):
        # notice key value must include #P or #p
        regex = re.compile('SCF Done')
        regex_1 = re.compile('Excitation Energies \[eV] at current iteration:')
        regex_2 = re.compile('Convergence achieved on expansion vectors|Convergence on energies')
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
                                    # get every states' energy/eV
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

    def grad(self, filename):
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
        return -np.array(data)

    def cico(self, time, filename):
        with open('traj_cicoe.log', 'a+') as f:
            regex = re.compile('Excited State\s*[0-9]*:')
            regex_1 = re.compile('SavETr')
            try:
                a = format(time, '>10d')
            except ValueError:
                a = format(time, '>10.2f') + '  fs'
            f.write("simulation time: t=" + a + '\n')
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

    def check(self, nstates):
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
                    if re.search(r'(?<=root=)%s' % nstates, line, re.IGNORECASE):
                        flag_r = True
                if not flag_p:
                    raise Exception("Please use detailed output(#P) ")
                if not flag_f:
                    raise Exception("Please use key value 'force'")
                if not flag_r:
                    raise Exception("Please check root=%s in gauss.gjf" % nstates)
        else:
            raise Exception("gauss.gjf is not found, please check *.gjf filename again")

    def r_wavefunction(self, filename):
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

    def d_wavefunction(self, filename):
        regex = re.compile('geom=allcheck', re.IGNORECASE)
        with open(filename, 'r') as f, open("%s.bak" % filename, 'w+') as f_new:
            for line in f:
                if not regex.search(line):
                    f_new.write(line)
                else:
                    f_new.write(re.sub(regex, '', line))
        os.remove(filename)
        os.rename("%s.bak" % filename, filename)

    def renew_calc_states(self, nstates, filename, filename_new=None, st=None, spin=None, charge=None, remove=None,
                          add=None):
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
                if flag_keyword:  # change the key word
                    if flag_td:
                        if nstates == 0:  # excited states -> ground states
                            line = re.sub(regex_td, '', line)
                        else:
                            line = re.sub(regex_root, str(nstates), line)
                    else:  # ground states -> the excited states
                        if re.search('#(p|P|s|S)?', line):
                            if nstates >= 1:
                                line = re.sub('\r?\n', '', line).strip()
                                WORD = re.sub(regex_root, str(nstates), _td_keyword)
                                line = line + ' ' + WORD + '\n'
                    if st:  # change the singlet/Triplets excited states for closed-shell systems
                        if re.search(regex_td, line):
                            if re.search(regex_st, line):
                                line = re.sub(regex_st, ST[st], line)
                            else:
                                LINE = re.search(
                                    "tda?=\(", line).group() + ST[st] + ","
                                line = re.sub("tda?=\(", LINE, line, re.IGNORECASE)
                    if remove:  # remove keywords
                        if re.search(remove, line, re.IGNORECASE):
                            line = re.sub(remove, "", line)
                    if add:  # add others keywords
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


def keyword(filename):
    regex = re.compile('(tda?=?\(.*?\)|tda)', re.IGNORECASE)
    key = []
    with open(filename) as f:
        for line in f:
            if regex.search(line):
                keyword = ' '.join(str(x) for x in regex.findall(line))
                key.append(keyword)
    tda = ' '.join(str(x) for x in key)
    return tda.strip()


_td_keyword = keyword('gauss.gjf')