import numpy as np
import matplotlib.pyplot as plt
import convergence
import time

res_dict = {
    0: "V_ra",
    1: "V_la",
    2: "V_rv",
    3: "V_lv",
    4: "Q_ra_rv",
    5: "P_pa",
    6: "Q_rv_pa",
    7: "Q_la_lv",
    8: "P_ao",
    9: "Q_lv_ao",
    10: "P_sys"
}

g_dict = {
    0: "P_ra",
    1: "P_rv",
    2: "P_la",
    3: "P_lv",
    4: "Q_sys_v",
    5: "Q_sys_a",
    6: "Q_pa",
    7: "Phi_T",
    8: "Phi_P",
    9: "Phi_M",
    10: "Phi_A"
}

def test_convergence(model, taus, ref_tau=1/16000):
    '''

    :param model: model object to run simulation
    :param taus: array of taus to test
    :param param: which parameter to test
    :param set_g: boolean if parameter is in g instead of u
    :param ref_tau: fine grid tau
    :return: array of errors relative to fine grid
    '''

    no_params = 11

    u_convergences = [[] for x in range(0,no_params)]
    g_convergences = [[] for x in range(0,no_params)]
    t_disc_ref, u_ref, g_ref = model.solve(model.cycleTime * ref_tau)
    u_ref_params = [u_ref[:, param] for param in range(0, no_params)]
    g_ref_params = [g_ref[:, param] for param in range(0, no_params)]
    for tau in taus:
        st = time.time()
        print(tau, "is starting!")
        t_disc, u, g = model.solve(model.cycleTime * tau)
        u_params = [u[:, param] for param in range(0, no_params)]
        g_params = [g[:, param] for param in range(0, no_params)]
        for param in range(0, no_params):
            u_convergences[param].append(convergence.convergence_norm(u_params[param], tau, u_ref_params[param], ref_tau))
            g_convergences[param].append(convergence.convergence_norm(g_params[param], tau, g_ref_params[param], ref_tau))
        et = time.time()
        print(tau, "is done! Calculations took ", et - st, " seconds.")
    # print(u_convergences, g_convergences)
    return u_convergences, g_convergences

def plot_convergence(model):

    taus = [1/500, 1/1000, 1/2000, 1/4000, 1/8000]
    # taus = [1/500, 1/1000]
    default_x_ticks = range(len(taus))

    u_conv, g_conv = test_convergence(model, taus, ref_tau=1/16000)

    for i in range(0, len(u_conv)):
        plt.plot(default_x_ticks, u_conv[i], 'o-', label=res_dict[i])

    plt.xticks(default_x_ticks, taus)
    plt.xlabel(r'$\tau$')
    plt.ylabel("Error")
    plt.yscale('log')
    plt.legend()
    plt.title("Convergence for ODE parameters")
    plt.savefig("convergence u.svg")
    plt.show()

    for i in range(0, len(g_conv)):
        plt.plot(default_x_ticks, g_conv[i], 'o-', label=g_dict[i])

    plt.xticks(default_x_ticks, taus)
    plt.xlabel(r'$\tau$')
    plt.ylabel("Error")
    plt.yscale('log')
    plt.legend()
    plt.title("Convergence for Algebraic parameters")
    plt.savefig("convergence g.svg")
    plt.show()



