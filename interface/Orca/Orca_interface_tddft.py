import os
import re

import numpy as np
from typing import List, Tuple, Union

from ..PublicFunction import PublicFunction

eV = 27.21138602
ang = 0.529177257507  # Angstrom/Å e-10


class Orca(PublicFunction):
    def grad(self, filename: str) -> np.ndarray:
        regex = re.compile("CARTESIAN GRADIENT")
        data: List = []
        flag_b: bool = False
        flag_g: bool = False
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

    def energy(self, filename: str, nstates: int) -> Tuple[List, float]:
        # notice, "TDDFT module" must include printlevel 3
        regex_scf = re.compile("Total Energy       :")
        regex_l1 = re.compile("lowest eigenvalues of")
        regex_l2 = re.compile("Lowest Energy  ")
        regex_e = re.compile("E\(tot\)  =")
        flag_l1: bool = False
        flag_l2: bool = False
        flag_scf: bool = False
        e_ground: float = 0.0
        all_excited: float = []
        e_excited: float = 0.0
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

    def cico(self, time: Union[float, int], filename: str) -> None:
        with open("traj_cicoe.log", "a+") as f, open(filename, 'r', encoding='UTF-8') as f_o:
            regex = re.compile("TD-DFT(/TDA)?\s*EXCITED STATES")
            regex_end = re.compile("E\(tot\)  =")
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
                        f.write("-" * 36 + "\n")
                        f.write(line)
                else:
                    f.write(line)
                    if regex_end.search(line):
                        f.write("\n")
                        break

    def check(self, nstates: int):
        if os.path.isfile("orca.inp"):
            flag_f = False
            flag_r = False
            flag_p = False 
            flag_w = False # %moinp "orca.gbw"
            flag_b = False # %base "orca"
            with open("orca.inp", 'r') as f:
                for line in f:
                    if re.search("EnGrad", line, re.IGNORECASE):
                        flag_f = True
                    if re.search("iroot\s+(?P<state>[1-9])", line, re.IGNORECASE):
                        state_str = re.search("iroot\s+(?P<state>[1-9])",line,re.IGNORECASE).group("state")
                        if str(nstates) != state_str :
                            if nstates == 0:
                                raise Exception("Please remove tddft module")
                            else:
                                raise Exception("Please check 'iroot %s' in orca.inp" % nstates)
                    if re.search("printlevel\s+3", line, re.IGNORECASE):
                        flag_p = True
                    if re.search(r"moinp \"orca.old.gbw\"", line, re.IGNORECASE):
                        flag_w = True 
                    if re.search(r"base\s+\"orca\"", line, re.IGNORECASE):
                        flag_b = True 
                if not flag_f:
                    raise Exception("Please use keyword 'EnGrad'")
                if not flag_w:
                    raise Exception("Please use keyword 'moread'")
                if not flag_b:
                    raise Exception("Please use keyword %base \"orca\"")
                if nstates >= 1:
                    if not flag_p:
                        raise Exception("Please use key value 'printlevel 3'")
        else:
            raise Exception("orca.inp is not found, please check *.inp filename again")

    def r_wavefunction(self, filename):
        pass

    def d_wavefunction(self, filename):
        pass

    def renew_calc_states(self, nstates: int, filename: str, filename_new: str = None,
                          st=None, spin=None, charge=None, remove=None, add=None):
        # exclamation(!)!CAM-B3LYP nousesym pal4  D3 RIJCOSX def2-SVP def2/J def2-SVP/C miniprint  TightSCF grid4 EnGrad
        # %pal nprocs 4 end
        # %TDDFT NROOTS 5
        # IROOT    1
        # printlevel 3
        # tda false
        # end
        flag_td: bool = False
        flag_keyword = False
        flag_td_end: bool = False
        flag_file: bool = False
        regex_td = re.compile("%TDDFT", re.IGNORECASE)
        count: int = 0
        td_n: int = 0
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
                            old_line = re.sub(
                                'iroot\s+[0-9]+', "iroot  %s" % nstates, _td_keyword, flags=re.IGNORECASE)
                            line = old_line + line
                f_new.write(line)

        if flag_file:
            os.remove(filename)
            os.rename(filename_new, filename)


def keyword(filename):
    regex = re.compile("%TDDFT", re.IGNORECASE)
    key = []
    flag_td = False
    # %TDDFT Nroot 5
    #        iroot 1
    #        tda false
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
    return "".join(x for x in key)


_td_keyword = keyword('orca.inp')