import os
import re

import numpy as np

from typing import List, Tuple, Union
from ..PublicFunction import PublicFunction

eV = 27.21138602
ang = 0.529177257507  # Angstrom/Å e-10


class Molpro(PublicFunction):
    def energy(self, filename: str, nstates: int) -> Tuple[List, float]:
        # !MCSCF STATE 2.1 Energy             -591.497026758027
        # the is singlets energy
        regex = re.compile("!MCSCF STATE \d+\.\d+ Energy")
        data = []
        with open(filename, 'r') as f:
            for line in f:
                if regex.search(line):
                    data.append(float(line.split()[-1]))
        return sorted(data), sorted(data)[nstates]

    def grad(self, filename: str) -> np.ndarray:
        regex = re.compile("SA-MC GRADIENT FOR STATE")
        flag: bool = False
        data: List = []
        with open(filename, 'r') as f:
            blank_count = 0
            for line in f:
                if not flag:
                    if regex.search(line):
                        flag = True
                else:
                    if not line.strip():
                        blank_count += 1
                    else:
                        if blank_count == 2:
                            data.append(list(map(float, line.split()[1:4])))
                        elif blank_count >= 2:
                            break
        return np.array(data)

    def cico(self, time: Union[float, int], filename: str) -> None:
        regex_ci = re.compile("CI vector")
        regex_end = re.compile("TOTAL ENERGIES")
        with open("traj_cicoe.log", "a+") as f, open(filename, 'r', encoding="UTF-8") as f_o:
            flag = False
            try:
                a = format(time, '>10d')
            except ValueError:
                a = format(time, '>10.2f') + '  fs'
            f.write("simulation time: t=" + a + "\n")
            for line in f_o:
                if not flag:
                    if regex_ci.search(line):
                        flag = True
                        f.write("-" * 36 + "\n")
                        f.write(line)
                else:
                    f.write(line)
                    if regex_end.search(line):
                        f.write("\n")
                        break

    def r_wavefunction(self, filename):
        pass

    def d_wavefunction(self, filename):
        pass

    def renew_calc_states(self, nstates: int, filename: str, filename_new: str = None, spin: int = 0, **kwargs):
        # CPMCSCF,GRAD,1.1,spin=0,accu=1.0d-7,record=5101.1 !gradient for state 1
        # CPMCSCF,GRAD,state,[SPIN=spin],[MS2=ms2],[ACCU=thresh],[RECORD=record]
        # Force;SAMC,record
        # accu: [+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eEdD][+-]?[0-9]+)?  => 1.0d-7
        # state: [0-9]\.[0-9]  =>2.1
        # spin: ([0-9]|[0-9]\.[0-9]) => 0  or 0.5 or 1
        # record: [0-9]*\.[0-9] =>  5101.1
        grad_1 = "CPMCSCF,GRAD,(?P<state>[0-9]\.[0-9])\s*,spin=(?P<spin>[0-9]|[0-9]\.[0-9]),"
        grad_2 = "accu=(?P<accu>[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eEdD][+-]?[0-9]+)?),"
        record_str = "record=(?P<record>[0-9]*\.[0-9])"
        regex_grad = re.compile(grad_1 + grad_2 + record_str, re.IGNORECASE)
        regex_force = re.compile("{Force;SAMC," + "(?P<record>[0-9]*\.[0-9])" + "}", re.IGNORECASE)

        def cpmscf(state: str = "1.1", s: str = '0', accu: str = "1.0d-7", record: str = "5101.1") -> str:
            v = "CPMCSCF,GRAD," + state + "," + "spin=" + s + ","
            v = v + "accu=" + accu + "," + "record=" + record
            return v

        def force(record: str = "5101.1") -> str:
            return "{Force;SAMC," + record + "}"

        regex_wf = re.compile("WF,(?P<elec>[0-9]*),(?P<sym>[0-9]),(?P<spin>[0-9]|[0-9]\.[0-9])", re.IGNORECASE)
        flag_file: bool = False
        sy: str = '1'
        sp: str = spin
        record_new: str = ''
        if not filename_new:
            flag_file = True
            filename_new = "%s.bak" % filename
        with open(filename, 'r') as f, open(filename_new, 'w+') as f_new:
            for line in f:
                # ignore the notes line (!xxxxx)
                if not re.match('^\s*!', line):
                    if regex_wf.search(line):
                        a = regex_wf.search(line)
                        sy = a.group("sym")
                        sp = a.group("spin")
                    if regex_grad.search(line):
                        n = str(nstates + 1) + "." + str(sy)
                        record_new = "510" + str(nstates + 1) + ".1"
                        line = cpmscf(state=n, s=sp, record=record_new) + "\n"
                    if regex_force.search(line):
                        line = force(record=record_new) + "\n"
                f_new.write(line)
        if flag_file:
            os.remove(filename)
            os.rename(filename_new, filename)

    def check(self, nstates: int, spin: int = 0):
        # if os.path.isfile("./MolproTemp/molpro.wfu"):
        #     print("The wavefunction file 'molpro.wfu' exist in ./MolproTemp and you should delete it")
        #     os.remove("./MolproTemp/molpro.wfu")
        if os.path.isfile("molpro.in"):
            # file,2,**.wfn
            record_str: str = ''
            flag_wfn: bool = False
            flag_grad: bool = False
            flag_force: bool = False
            with open("molpro.in", 'r') as f:
                for line in f:
                    if not re.match('^\s*!', line):
                        # ignore black notes
                        if re.search("file,2,.*\.wfu", line, re.IGNORECASE):
                            flag_wfn = True
                        if re.search("^CPMCSCF,GRAD,%s\.1,spin=%s" % (str(nstates + 1),
                                                                      str(spin)), line, re.IGNORECASE):
                            record_str = re.search("record=(?P<record>[0-9]*\.[0-9])", line).group("record")
                            # record default is 510N.1 ,N is current states
                            if not (record_str == '510' + str(nstates + 1) + '.1'):
                                raise Exception("Record(record=%s) maybe is error,Please check 'molpro.in' again"
                                                % record_str)
                            flag_grad = True
                        if re.search("{Force;SAMC,%s}" % record_str, line, re.IGNORECASE):
                            flag_force = True
                if not flag_wfn:
                    raise Exception("Please read wfu file in 'file,2,molpro.wfu' ")
                if not flag_grad:
                    raise Exception(
                        "Please check ’CPMCSCF,GRAD,%s.1,spin=%s....’ in molpro.in" % (str(nstates + 1), str(spin)))
                if not flag_force:
                    raise Exception("Please check '{Force;SAMC,%s}' in molpro.in" % record_str)
        else:
            raise Exception("molpro.in is not found, please check *.in filename again")