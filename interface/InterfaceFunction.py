from utils.geom import replace_coord
from .PublicFunction import PublicFunction
from .Gaussian import Gaussian
from .Orca import Orca
from .prepare import (soft, dyn_states,
                      output_suffix, input_suffix, input_prefix)

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

    def grad(self, filename=None):
        if not filename:
            filename = self.output
        return self.result.grad(filename)

    def energy(self, filename=None):
        if not filename:
            filename = self.output
        return self.result.energy(filename)

    def cico(self, time, filename=None):
        if not filename:
            filename = self.output
        return self.result.cico(time, filename)

    def replace(self, new_coord, filename=None, filename_new=None):
        if not filename:
            filename = self.input
        return replace_coord(new_coord, filename, filename_new)

    def renew_states(self, nstates, filename=None, filename_new=None, **kwargs):
        if not filename:
            filename = self.input
        return self.result.renew_calc_states(nstates, filename, filename_new, **kwargs)

    def read(self, filename=None):
        if not filename:
            filename = self.input
        return self.result.r_wavefunction(filename)

    def delete(self, filename=None):
        if not filename:
            filename = self.input
        return self.result.d_wavefunction(filename)

    def chk(self, nstates):
        return self.result.check(nstates)


s = eval(soft + "()")
r = Interface(s)


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
    replace new coordinate in inputfile
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
    Check whether the initial conditions are reasonable
    and whether the inputfile( "gauss.gjf, molpro.in, orca.inp") keywords is wrong
    Args:
        nstates: the current dynamic states

    Returns:

    """

    return r.chk(nstates)


check_initial(dyn_states)


def renew_calc_states(nstates, filename=None, filename_new=None,
                      st=None, spin=None, charge=None, remove=None, add=None):
    """
    Calculate the energy of different states(S1->S0,S0-S1, S->ST(singlets,triplets)
    change some keywords(spin, charges, add and delete keyword)
    Args:
        nstates:  the new dynamic states
        filename:  default is "gauss.gjf, orca.inp, molpro.in"
        filename_new: default is "None"
        st: (singlets,triplets) is only calculate energy
        spin: singles, triplets
        charge: 0 1
        remove: remove the old keyword in last line
        add: add some new keyword in last line
    Returns:
        None
    """
    return r.renew_states(nstates, filename, filename_new,
                          st=None, spin=None, charge=None, remove=None, add=None)