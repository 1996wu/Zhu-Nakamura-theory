from typing import List, Tuple, Union
import numpy as np

from utils.geom import replace_coord
from .PublicFunction import PublicFunction
from .prepare import (soft, dyn_states,
                      output_suffix, input_suffix, input_prefix)

if soft == "Gaussian":
    from .Gaussian import Gaussian
elif soft == "Orca":
    from .Orca import Orca
elif soft == "Molpro":
    from .Molpro import Molpro

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

    def energy(self, filename=None, nstates=None):
        if not filename:
            filename = self.output
        return self.result.energy(filename, nstates)

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


def get_energy(filename: str = None, nstates: int = None) -> Tuple[np.ndarray, float]:
    """

    Args:
        filename: The output file of quantum chemistry software: "gauss.log", "orca.out", "molpro.out":
        nstates: the energy of the current states(This must be pointed in CASSCF level with Molpro)
    Returns:
        tuple,The excited energy list [S0, S1, S2, S3] and  the current total energy()
    """
    return r.energy(filename, nstates)


def get_grad_matrix(filename: str = None) -> np.ndarray:
    """

    Args:
        filename: The output file of quantum chemistry software: "gauss.log", "orca.out" ...

    Returns:
       The grad matrix(3* Natom)
    """
    return r.grad(filename)


def print_traj_cicoe(time: Union[float, int], filename: str = None):
    """

    Args:
        time: The current dynamics time(default float) or The run time(int)
        filename:  The output file of quantum chemistry software: "gauss.log", "orca.out" ...

    Returns: None
    """
    return r.cico(time, filename)


def replace_coordinate(new_coord: np.ndarray, filename: str = None, filename_new: str = None):
    """
    replace new coordinate in inputfile
    Args:
        new_coord: the geom coordinate(a.u.)(3 * N)
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


def check_initial(nstates: int) -> None:
    """
    Check whether the initial conditions are reasonable
    and whether the inputfile( "gauss.gjf, molpro.in, orca.inp") keywords is wrong
    Args:
        nstates: the current dynamic states

    Returns:
            None
    """
    return r.chk(nstates)


check_initial(dyn_states)


def renew_calc_states(nstates: int, filename: str = None, filename_new: str = None,
                      **kwargs):
    """
    Calculate the energy of different states(S1->S0,S0-S1, S->ST(singlets,triplets)
    change some keywords(spin, charges, add and delete keyword)
    Args:
        nstates:  the new dynamic states
        filename:  default is "gauss.gjf, orca.inp, molpro.in"
        filename_new: default is "None"
        Kwargs:
            st: (singlets,triplets) is only calculate energy
            spin: singles, triplets
            charge: 0 1
            remove: remove the old keyword in last line
            add: add some new keyword in last line
    Returns:
        None
    """
    return r.renew_states(nstates, filename, filename_new,**kwargs)