#!/usr/bin/env python

"""
This is interface with GAUSSIAN, ORCA
include: 
software_running,read_wavefunction, 
delete_wavefunction, check_initial,
get_grad_matrix, get_energy, 
replace_coordinate, renew_calc_states,
print_traj_cicoe
Author: zibo wu <zbwu1996@gmail.com>
"""

from .InterfaceFunction import *
from .RunSoftware import software_running

__all__ = ["software_running", "read_wavefunction", "delete_wavefunction", "check_initial",
           "get_grad_matrix", "get_energy", "replace_coordinate", "renew_calc_states", "print_traj_cicoe",
           "eV", "ang", "fs", "amu", "pi", "velo"]

eV = 27.21138602
ang = 0.529177257507  # Angstrom/Å e-10
fs = 0.0241888439241  # 1 a.u. = 2.4188*10e-17
amu = 1822.88853006  # Relative atomic mass
pi = 3.141592653589793
velo = 2.18769126364 * 10e6  # 1 a.u. = 2.187*10e4  m/s