import copy
import re
import os 
import numpy as np

from interface.prepare import flag_rot, natom, element_mass, element

__all__ = ["critical_value", "coord_corrections", "mom_corrections",
           "calculate_kinetic", "print_matrix", "replace_coord"]
ang = 0.529177257507


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


def replace_coord(new_coordinate, filename, filename_new=None):
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
                new_coord = ''.join(format(i * ang, '>18.10f')
                                    for i in new_coordinate[count][:3])
                f_new.write(re.sub(regex_1, new_coord, line))
                count += 1
        if count != len(new_coordinate):
            raise Exception("the atom of file(%s, %s) is not equal of new_coordinate(%s)"
                                %(filename, str(count), str(len(new_coordinate))))
    if flag_file:
        os.remove(filename)
        os.rename(filename_new, filename)


def calculate_kinetic(mom):
    kine = 0.0
    for i in range(natom):
        kine += 0.5 * np.sum(mom[i] ** 2) / element_mass[i]
    return kine


def print_matrix(value):
    for i in range(natom):
        ele = format(element[i], '<5s')
        xyz = "".join(format(x, '>18.10f') for x in value[i])
        print(ele + xyz, flush=True)