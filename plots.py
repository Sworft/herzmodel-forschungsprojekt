import matplotlib.pyplot as plt

# mmHg to Barye Conversion Factor
mmHgToBarye = 1333.22
baryeTommHg = 1.0 / 1333.22


def overview_plots(model, tau, cycleTime = 1.07, show_volume=1, show_flow=1, show_pressure=1, show_valve=1):

    t_disc, u, g = model.solve(cycleTime * tau)

    for i in range(len(u)):  # convert Barye pressure to mmHg
        u[i][5] *= baryeTommHg
        u[i][8] *= baryeTommHg
        u[i][10] *= baryeTommHg
        g[i][0] *= baryeTommHg
        g[i][1] *= baryeTommHg
        g[i][2] *= baryeTommHg
        g[i][3] *= baryeTommHg


    if show_volume:
        fig = plt.figure()
        ax = fig.add_axes([0.2, 0.25, 0.7, 0.6])
        ax.set_title('Volume', fontsize=25)
        ax.set_xlabel('Time', fontsize=25)
        ax.set_ylabel('[mL]', fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=20)

        line, = ax.plot(t_disc, u[:, 0], label="V_ra")
        line, = ax.plot(t_disc, u[:, 1], label="V_la")
        line, = ax.plot(t_disc, u[:, 2], label="V_rv")
        line, = ax.plot(t_disc, u[:, 3], label="V_lv")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35),
                  fancybox=True, shadow=True, ncol=5, prop={'size': 13})

        plt.show()

    if show_flow:
        fig = plt.figure()
        ax = fig.add_axes([0.2, 0.25, 0.7, 0.6])
        ax.set_title('Flow rate', fontsize=25)
        ax.set_xlabel('Time', fontsize=25)
        ax.set_ylabel('[mL/s]', fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=20)

        # line, = ax.plot(t_disc, u[:, 4], label="Q_ra_rv")
        # line, = ax.plot(t_disc, u[:, 6], label="Q_rv_pa")
        # line, = ax.plot(t_disc, u[:, 7], label="Q_la_lv")
        # line, = ax.plot(t_disc, u[:, 9], label="Q_lv_ao")
        line, = ax.plot(t_disc, g[:, 4], label="Q_sys_v")
        line, = ax.plot(t_disc, g[:, 5], label="Q_sys_a")
        line, = ax.plot(t_disc, g[:, 6], label="Q_pa")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35),
                  fancybox=True, shadow=True, ncol=5, prop={'size': 13})

        plt.show()

    if show_pressure:

        fig = plt.figure()
        ax = fig.add_axes([0.2, 0.25, 0.7, 0.6])
        ax.set_title('Pressure', fontsize=25)
        ax.set_xlabel('Time', fontsize=25)
        ax.set_ylabel('[mmHG]', fontsize=25)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=20)

        # line, = ax.plot(t_disc, u[:, 5], label="P_pa")
        # line, = ax.plot(t_disc, u[:, 8], label="P_ao")
        # line, = ax.plot(t_disc, u[:, 10], label="P_sys")
        line, = ax.plot(t_disc, g[:, 0], label="P_ra")
        line, = ax.plot(t_disc, g[:, 1], label="P_rv")
        line, = ax.plot(t_disc, g[:, 2], label="P_la")
        line, = ax.plot(t_disc, g[:, 3], label="P_lv")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35),
                  fancybox=True, shadow=True, ncol=5, prop={'size': 13})

        plt.show()

    if show_valve:

        fig = plt.figure()
        ax = fig.add_axes([0.2, 0.25, 0.7, 0.6])
        ax.set_title('Valve state', fontsize=25)
        ax.set_xlabel('Time', fontsize=25)
        ax.set_ylabel("", fontsize=25)
        plt.yticks([0,1], ["Closed","Open"])
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=20)

        line, = ax.plot(t_disc, g[:, 7], label="Phi_T")
        line, = ax.plot(t_disc, g[:, 8], label="Phi_P")
        line, = ax.plot(t_disc, g[:, 9], label="Phi_M")
        line, = ax.plot(t_disc, g[:, 10], label="Phi_A")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35),
                  fancybox=True, shadow=True, ncol=5, prop={'size': 13})
        plt.show()