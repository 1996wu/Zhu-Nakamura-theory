import numpy as np
import copy

from interface.prepare import flag_rot, natom, element_mass

__all__ = ["critical_value", "coord_mom_rot", "coord_corrections", "mom_corrections",
           "calculate_kinetic"]


def critical_value(*args):
    def __dihedral_angle(p1, p2, p3, p4):
        q1 = np.subtract(p2, p1)
        q2 = np.subtract(p3, p2)
        q3 = np.subtract(p4, p3)
        # Calculate cross vectors
        q1_x_q2 = np.cross(q1, q2)
        q2_x_q3 = np.cross(q2, q3)
        # Calculate normal vectors
        n1 = q1_x_q2 / np.linalg.norm(q1_x_q2)
        n2 = q2_x_q3 / np.linalg.norm(q2_x_q3)
        # Calculate unit vectors
        u1 = n2
        u3 = q2 / np.linalg.norm(q2)
        u2 = np.cross(u3, u1)
        # Calculate cosine and sine
        cos_theta = np.dot(n1, u1)
        sin_theta = np.dot(n1, u2)
        theta = -np.arctan2(sin_theta, cos_theta)
        theta_deg = np.degrees(theta)
        return theta_deg

    def __bond_angle(p1, p2, p3):
        q1 = np.subtract(p2, p1)
        q2 = np.subtract(p3, p2)
        cos_theta = np.dot(q1, q2) / (np.linalg.norm(q1) * np.linalg.norm(q2))
        theta = np.arccos(cos_theta)
        theta_deg = 180.00 - np.degrees(theta)
        return theta_deg

    def __bond_length(p1, p2):
        length = np.linalg.norm(np.subtract(p2, p1))
        return length

    N = len(args)
    if N == 2:
        return __bond_length(*args)
    elif N == 3:
        return __bond_angle(*args)
    elif N == 4:
        return __dihedral_angle(*args)
    else:
        raise Exception("Too many parameters")


def coord_mom_rot(coord1, coord2):
    """
    Kabsch algorithm
    url=https://en.wikipedia.org/wiki/Kabsch_algorithm
    """
    if not flag_rot:
        return np.eye(3)
    # natom = len(coord1)
    coord_1 = copy.deepcopy(coord1)
    coord_2 = copy.deepcopy(coord2)
    coord_1 -= np.mean(coord_1, axis=0)
    coord_2 -= np.mean(coord_2, axis=0)
    cov_matrix = coord_1.T @ coord_2
    U, sigma, VT = np.linalg.svd(cov_matrix)
    d = np.sign(np.linalg.det(U @ VT))
    diag = np.identity(3)
    diag[2][2] = d
    rot_mat = VT.T @ diag @ U.T
    # rot_coord = coord_2 @ rot_matrix
    # RMSD = np.sqrt(np.sum((coord_1 - rot_coord) ** 2) / natom)
    # 消除平移 https://zhuanlan.zhihu.com/p/111322916
    return rot_mat


def coord_corrections(coord0, geom, mass):
    """
    Modifies the coordinate to center of mass
    geom: N * 3;mass : 1 * N
    """
    total_mass = np.sum(np.array(mass))
    mass_N = np.array(mass).reshape((len(mass), 1)) / total_mass * len(mass)
    coord_rot_matrix = coord_mom_rot(coord0 * mass_N, geom * mass_N)
    new_geom = geom @ coord_rot_matrix
    return new_geom, coord_rot_matrix


def mom_corrections(mom1, mom2):
    """
    mom: N * 3
    """
    mom_rot_matrix = coord_mom_rot(mom1, mom2)
    new_mom = mom2 @ mom_rot_matrix
    return new_mom


def calculate_kinetic(mom):
    kine = 0.0
    for i in range(natom):
        kine += 0.5 * np.sum(mom[i] ** 2) / element_mass[i]
    return kine