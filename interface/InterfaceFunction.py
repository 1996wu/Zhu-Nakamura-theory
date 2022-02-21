import re
import os
import numpy as np

from .Orca import Orca
from .Gaussian import Gaussian
from .PublicFunction import PublicFunction
from .prepare import soft, dyn_states
from .prepare import output_suffix, input_suffix, input_prefix

ang = 0.529177210903
Soft_class = {"Gaussian": "Gaussian()", "Orca": "Orca()", "Molpro": "Molpro()",
              "Bdf": "Bdf()"}

__all__ = ["read_wavefunction", "delete_wavefunction", "check_initial",
           "get_grad_matrix", "get_energy", "replace_coordinate", "renew_calc_states", "print_traj_cicoe"]


class Interface:

    def __init__(self, c: PublicFunction):
        self.soft = soft
        self.output = input_prefix[soft] + "." + output_suffix[soft]
        self.input = input_prefix[soft] + "." + input_suffix[soft]
        self.result = c

    @staticmethod
    def __replace(new_coordinate, filename, filename_new=None):
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

    def grad(self, filename=None):
        if filename:
            self.output = filename
        return self.result.grad(self.output)

    def energy(self, filename=None):
        if filename:
            self.output = filename
        return self.result.energy(self.output)

    def cico(self, time, filename=None):
        if filename:
            self.output = filename
        return self.result.cico(time, self.output)

    def replace(self, new_coord, filename=None, filename_new=None):
        if filename:
            self.input = filename
        return self.__replace(new_coord, self.input, filename_new)

    def renew_states(self, nstates, filename=None, filename_new=None, **kwargs):
        if filename:
            self.input = filename
        return self.result.renew_calc_states(nstates, self.input, filename_new, **kwargs)

    def get_keyword(self, filename=None):
        if filename:
            self.input = filename
        return self.result.keyword(self.input)

    def read(self, filename=None):
        if filename:
            self.input = filename
        return self.result.r_wavefunction(self.input)

    def delete(self, filename=None):
        if filename:
            self.input = filename
        return self.result.d_wavefunction(self.input)

    def chk(self, nstates):
        return self.result.check(nstates)


s = eval(soft + "()")
r = Interface(s)
word = r.get_keyword(filename=None)


def get_energy(filename=None):
    """

    Args:
        filename: The output file of quantum chemistry software: "gauss.log", "orca.out" ..:

    Returns:
        The excited energy list [S0, S1, S2, S3] and  the current total energy(float)
    """
    return r.energy(filename)


def get_grad_matrix(filename=None):
    """

    Args:
        filename: The output file of quantum chemistry software: "gauss.log", "orca.out" ...

    Returns:
       The grad matrix(3* Natom)
    """
    return r.grad(filename)


def print_traj_cicoe(time, filename=None):
    """

    Args:
        time: The current dynamics time(default float) or The run time(int)
        filename:  The output file of quantum chemistry software: "gauss.log", "orca.out" ...

    Returns: None
    """
    return r.cico(time, filename)


def replace_coordinate(new_coord, filename=None, filename_new=None):
    """

    Args:
        new_coord: the geom coordinate(a.u.)
        filename:  default is (gauss.gjf, orca.inp, molpro.in)
        filename_new: Whether to save as a new file, default is "None"

    Returns:
            None
    """
    return r.replace(new_coord, filename, filename_new)


def read_wavefunction(filename=None):
    return r.read(filename)


def delete_wavefunction(filename=None):
    return r.delete(filename)


def check_initial(nstates):
    """

    Args:
        nstates: the current dynamic states

    Returns:

    """

    return r.chk(nstates)


check_initial(dyn_states)


def renew_calc_states(nstates, filename=None, filename_new=None, **kwargs):
    """

    Args:
        nstates:  the new dynamic states
        filename:  default is "gauss.gjf, orca.inp, molpro.in"
        filename_new: default is "None"
        **kwargs: : add=None(add the new keywords), remove=None(remove the old keyword),spin

    Returns:

    """
    return r.renew_states(nstates, filename, filename_new, **kwargs)