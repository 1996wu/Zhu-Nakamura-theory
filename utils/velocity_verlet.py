import numpy as np

from interface.prepare import SI_step_time, fs, natom
from interface.prepare import element_mass

__all__ = ["update_position_matrix", "update_momentum_matrix"]


def update_position_matrix(position_matrix, momentum_matrix, grad_matrix,
                           mass=element_mass, delta=SI_step_time):
    """
    x_n+1 = x_n +  p_n /m_i * Δt  + 1/2mi  * F_n * Δt^2
    """
    step_time = delta / fs
    position_matrix_back = []
    for i in range(natom):
        data = position_matrix[i] + (momentum_matrix[i] / mass[i]) * step_time - 0.5 * \
               (grad_matrix[i] / mass[i]) * step_time ** 2
        position_matrix_back.append(data)
    return np.array(position_matrix_back)


# notice atom force symbol positive + or negative  -

def update_momentum_matrix(momentum_matrix, grad_matrix, grad_matrix_back, delta=SI_step_time):
    """
    p_n+1 = p_n + 0.5 * (F_n+1 +F_n) * Δt
    """
    step_time = delta / fs
    velocity_matrix_back = []
    for i in range(natom):
        data = momentum_matrix[i] - 0.5 * (grad_matrix[i] + grad_matrix_back[i]) * step_time
        velocity_matrix_back.append(data)
    return np.array(velocity_matrix_back)