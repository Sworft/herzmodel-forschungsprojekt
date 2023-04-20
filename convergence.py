import numpy as np

def prolongate(u, tau, ref_tau = 1/16000):
    """

    :param u: coarse grid results to be prolongated
    :param tau: coarse grid tau
    :param ref_tau: fine grid reference tau
    :return: prolongated result vector with linearly interpolated values as numpy array
    """
    step_size = int(tau/ref_tau)
    u_p = []
    for i in range(0, len(u)-1):
        m = u[i+1] - u[i]
        for j in range(0, step_size):
            u_p.append(u[i] + m * j/step_size)
    u_p.append(u[-1])
    return np.array(u_p)

def trapez(u,t,ti):
    """

    :param u: function/vector to be integrated
    :param t: timestep i
    :param ti: timestep i+1
    :return: approximation of integral using trapez method as float
    """
    return (ti - t) * (1/2) * (u[t]**2 + u[ti]**2)

def l2(u):
    """

    :param u: function/vector to be normed
    :return: L2 norm, approximated using trapez method as float
    """
    norm = 0
    for i in range(0, len(u)-1):
        norm += trapez(u, i, i+1)
    return norm

def convergence_norm(u, tau, u_ref, ref_tau=1/16000):
    """

    :param ref_tau: fine grid tau
    :param u: coarse grid results
    :param tau: coarse grid tau
    :param u_ref: fine grid reference
    :return: convergence as L2 norm of the difference between coarse and fine grid solutions normalized to fine grid as float
    """
    u_p = prolongate(u, tau, ref_tau)
    u_dif = u_ref - u_p
    norm_dif = l2(u_dif)
    norm_ref = l2(u_ref)
    return norm_dif / norm_ref

    # return np.linalg.norm(u_dif)/np.linalg.norm(u_ref)