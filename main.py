import numpy as np
import matplotlib.pyplot as plt
import analysis
import plots

# mmHg to Barye Conversion Factor
mmHgToBarye = 1333.22
baryeTommHg = 1.0 / 1333.22

# labels

def getMean(t,y):
  '''
  Get the mean of a signal in time
  Computing the integral and dividing by the total time
  '''
  return np.trapz(y,t)/(t[-1]-t[0])

def explicit_euler(f, u0, g0, param, t_disc):
    '''
    Forward Euler integration method
    :param f: system of differential equations
    :param u0: ODE outputs
    :param g0: Algebraic outputs
    :param param: Model parameters
    :param t_disc: Time discretization
    :return: Arrays containing each time step (ODE/Alg)
    '''

    # Neues array mit u0 als erstem Element.
    u = np.zeros((len(t_disc), len(u0)))
    g = np.zeros((len(t_disc), len(g0)))

    #Initialisierung \n",
    u[0] = u0 #die Unbekannten
    g[0] = g0 #die Algebraischen 'Unbekannten'

    #Zeitschritte durch für-Schleif:
    for i in range(1,len(t_disc)):
        t_last = t_disc[i-1]
        t = t_disc[i]
        tau = t - t_last

        # der letzte berechnete Wert von u
        u_last = u[i-1]
        g_last = g[i-1]

        # Ausführung des Zeitschritts
        du_dt, g_new = f(t_last, u_last, param, tau)
        u_new = u_last + tau * du_dt

        # Speichern des neuen Werts
        u[i] = u_new
        g[i] = g_new

    return u, g

def DAE(t, u, param, tau):

    # Algebraische Werte aus g

    # Herzfrequenz

    HR = param[0]

    # Atrial and ventricular duration and shift

    tsas = param[1]
    tpws = param[2]
    tsvs = param[3]

    # Atrial model parameters

    K_pas_ra_1 = param[4]
    K_pas_ra_2 = param[5]
    Emax_ra = param[6]
    Vra0 = param[7]

    K_pas_la_1 = param[8]
    K_pas_la_2 = param[9]
    Emax_la = param[10]
    Vla0 = param[11]

    #Ventricular Model Parameters

    K_pas_rv_1 = param[12]
    K_pas_rv_2 = param[13]
    Emax_rv = param[14]
    Vrv0 = param[15]
    K_pas_lv_1 = param[16]
    K_pas_lv_2 = param[17]
    Emax_lv = param[18]
    Vlv0 = param[19]

    #Atrial and Ventricular Inductances and Resistances

    L_ra = param[20]
    R_ra = param[21]
    L_rv = param[22]
    R_rv = param[23]
    L_la = param[24]
    R_la = param[25]
    L_lv = param[26]
    R_lv = param[27]

    # Aortic Arch

    C_ao = param[28]

    # Pulmonary Resistance and Capacitance

    C_pa = param[29]
    R_pa = param[30]

    # Systemic Resistance and Capacitance

    C_sys = param[31]
    R_sys_a = param[32]
    R_sys_v = param[33]





    # Letzte Werte von u

    #Volumen

    V_ra = u[0]
    V_la = u[1]
    V_rv = u[2]
    V_lv = u[3]

    #Flüsse

    Q_ra_rv = u[4]
    Q_rv_pa = u[6]
    Q_la_lv = u[7]
    Q_lv_ao = u[9]

    #Drücke

    P_pa = u[5]
    P_ao = u[8]
    P_sys = u[10]

    #Errechnete Algebraische Werte

    # Zeitvariablen

    tc = 60.0 / HR
    tsa = tc * tsas
    tsv = tc * tsvs
    tpw = tc / tpws
    tmv = t - np.floor(t/tc)*tc
    tma = (t+tsa-tpw) - np.floor((t+tsa-tpw)/tc)*tc
    fAV = (1 - np.cos(2 * np.pi * tmv / tsv)) / 2 if tmv < tsv else 0
    fAA = (1 - np.cos(2 * np.pi * tma / tsa)) / 2 if tma < tsa else 0


    #Drücke

    # ATRIA
    # Compute exponential atrial passive pressures curves

    P_pas_ra = K_pas_ra_1 * (np.exp((V_ra - Vra0) * K_pas_ra_2) - 1)
    P_pas_la = K_pas_la_1 * (np.exp((V_la - Vla0) * K_pas_la_2) - 1)

   # Compute linear atrial active pressure curves

    P_act_ra = Emax_ra * (V_ra - Vra0)
    P_act_la = Emax_la * (V_la - Vla0)

    # Blend with activation function

    P_ra = (P_pas_ra + fAA * (P_act_ra - P_pas_ra)) * mmHgToBarye
    P_la = (P_pas_la + fAA * (P_act_la - P_pas_la)) * mmHgToBarye

    # VENTRICLES
    # Passive Curve - Exponential

    P_pas_rv = K_pas_rv_1 * (np.exp((V_rv - Vrv0) * K_pas_rv_2) - 1.0)
    P_pas_lv = K_pas_lv_1 * (np.exp((V_lv - Vlv0) * K_pas_lv_2) - 1.0)

    # Active Curve - Linear

    P_act_rv = Emax_rv * (V_rv - Vrv0)
    P_act_lv = Emax_lv * (V_lv - Vlv0)

    # Use Activation to blend between active and passive Curves

    P_rv = (P_pas_rv + fAV * (P_act_rv - P_pas_rv)) * mmHgToBarye
    P_lv = (P_pas_lv + fAV * (P_act_lv - P_pas_lv)) * mmHgToBarye

    #Klappen

    Phi_T = 1 if P_ra > P_rv and Q_ra_rv >= 0 else 0
    Phi_P = 1 if P_rv > P_pa and Q_rv_pa >= 0 else 0
    Phi_M = 1 if P_la > P_lv and Q_la_lv >= 0 else 0
    Phi_A = 1 if P_lv > P_ao and Q_lv_ao >= 0 else 0

    # Flüsse
    Q_sys_a = (P_ao - P_sys) / R_sys_a
    Q_sys_v = (P_sys - P_ra) / R_sys_v
    Q_pa = (P_pa - P_la) / R_pa

    #Calculate ODE Rhs'

    #Volumen
    dV_ra = Q_sys_v - Q_ra_rv * Phi_T #V_ra: Volume of blood in right atrium [ml]
    dV_la = Q_pa - Q_la_lv * Phi_M #V_la: Volume of blood in left atrium [ml]
    dV_rv = Q_ra_rv * Phi_T - Q_rv_pa * Phi_P #V_rv: Volume of blood in right ventricle [ml]
    dV_lv = Q_la_lv * Phi_M - Q_lv_ao * Phi_A #V_lv: Volume of blood in left ventricle [ml]

    #Fluss

    # !!!!!!!! HIER IM GEGENSATZ ZUM PAPER REPARIERT !!!!!!!!!!!!
    # Im Paper wird 0 verwendet, falls die Klappen zu sind. Dies führt jedoch dazu, dass die Größen immer weiter anwachsen!
    # Ich habe es mit -Q_ra_rv/tau ersetzt, somit werden die Flüsse komplett auf 0 gesetzt, sobald die Klappen schließen.
    # Dies ist kongruenter mit der physiologischen Realität (wenn auch sehr vereinfacht) und mit dem Paper selbst (!),
    # wo dies genauso beschrieben wird, aber nicht implementiert oder in den Gleichungen aufgeführt.
    dQ_ra_rv = (P_ra - P_rv - R_ra * Q_ra_rv) / L_ra if Phi_T else -Q_ra_rv/tau #Q_ra,rv: Flow going from right atrium to right ventricle [ml/s]
    dQ_rv_pa = (P_rv - P_pa - R_rv * Q_rv_pa) / L_rv if Phi_P else -Q_rv_pa/tau #Q_rv,pa: Flow going from right ventricle to pulmonary artery [ml/s]
    dQ_la_lv = (P_la - P_lv - R_la * Q_la_lv) / L_la if Phi_M else -Q_la_lv/tau #Q_la,lv: Flow going from left atrium to left ventricle [ml/s]
    dQ_lv_ao = (P_lv - P_ao - R_lv * Q_lv_ao) / L_lv if Phi_A else -Q_lv_ao/tau #Q_lv,ao: Flow going from left ventricle to aortic arc [ml/s]

    # dQ_ra_rv = (P_ra - P_rv - R_ra * Q_ra_rv) / L_ra if Phi_T else 0 #Q_ra,rv: Flow going from right atrium to right ventricle [ml/s]
    # dQ_rv_pa = (P_rv - P_pa - R_rv * Q_rv_pa) / L_rv if Phi_P else 0 #Q_rv,pa: Flow going from right ventricle to pulmonary artery [ml/s]
    # dQ_la_lv = (P_la - P_lv - R_la * Q_la_lv) / L_la if Phi_M else 0 #Q_la,lv: Flow going from left atrium to left ventricle [ml/s]
    # dQ_lv_ao = (P_lv - P_ao - R_lv * Q_lv_ao) / L_lv if Phi_A else 0 #Q_lv,ao: Flow going from left ventricle to aortic arc [ml/s]


    #Druck
    dP_pa = (Q_rv_pa * Phi_P - Q_pa) / C_pa #P_pa: Pressure in pulmonary artery [Barye]
    dP_ao = (Q_lv_ao - Q_sys_a) / C_ao #P_ao: Pressure in aortic arc [Barye]
    dP_sys = (Q_sys_a - Q_sys_v) / C_sys #P_sys: Pressure in system [Barye]

    #Assign outputs

    # ODE_RHS

    res = np.zeros(len(u))


    #Volumen
    res[0] = dV_ra #V_ra: Volume of blood in right atrium [ml]
    res[1] = dV_la #V_la: Volume of blood in left atrium [ml]
    res[2] = dV_rv #V_rv: Volume of blood in right ventricle [ml]
    res[3] = dV_lv #V_lv: Volume of blood in left ventricle [ml]

    #Fluss
    res[4] = dQ_ra_rv #Q_ra,rv: Flow going from right atrium to right ventricle [ml/s]
    res[6] = dQ_rv_pa #Q_rv,pa: Flow going from right ventricle to pulmonary artery [ml/s]
    res[7] = dQ_la_lv #Q_la,lv: Flow going from left atrium to left ventricle [ml/s]
    res[9] = dQ_lv_ao #Q_lv,ao: Flow going from left ventricle to aortic arc [ml/s]

    #Druck
    res[5] = dP_pa #P_pa: Pressure in pulmonary artery [mmHg]
    res[8] = dP_ao #P_ao: Pressure in aortic arc [mmHg]
    res[10] = dP_sys #P_sys: Pressure in system [mmHg]






    # Algebraische Werte

    g = np.zeros(11)

    #Drücke
    g[0] = P_ra
    g[1] = P_rv
    g[2] = P_la
    g[3] = P_lv

    #Flüsse
    g[4] = Q_sys_v
    g[5] = Q_sys_a
    g[6] = Q_pa

    #Klappen
    g[7] = Phi_T
    g[8] = Phi_P
    g[9] = Phi_M
    g[10] = Phi_A

    return res, g

class HeartModel:

    def __init__(self, cycleTime, totalCycles, forcing=None, debugMode=False):
        # Init parameters
        self.u = None
        self.g = None
        self.cycleTime = cycleTime
        self.totalCycles = totalCycles
        self.tau = None
        self.post = None
        self.t_disc = None

        # NOTE: CGS Units: Pressures in Barye, Flowrates in mL/s
        # Default Initial Conditions
        self.defIC = np.array([  # Initial Values
            0.0,  # V_ra
            0.0,  # V_la
            0.0,  # V_rv
            0.0,  # V_lv
            0.0,  # Q_ra_rv
            70.0 * mmHgToBarye,  # P_pa
            0.0,  # Q_rv_pa
            0.0,  # Q_la_lv
            100.0 * mmHgToBarye,  # P_ao
            0.0,  # Q_lv_ao
            50.0 * mmHgToBarye])  # P_syss

        # NOTE: CGS Units: Pressures in Barye, Flowrates in mL/s
        self.defParam = np.array([  # Heart Cycle Parameters
            78.0,  # HR - Heart Rate (beats per minute)
            # Atrial and ventricular activation duration
            0.2,  # tsas - Atrial relative activation duration
            9.5,  # tpws - Atrial relative activation time shift
            0.4,  # tsvs - Ventricular relative activation duration
            # Atrial model parameters
            5.0,  # K_pas_ra_1 - Atrial passive curve slope, right atrium
            0.006,  # K_pas_ra_2 - Atrial passive curve exponent factor, right atrium
            0.1,  # Emax_ra - Atrial active curve slope, right atrium
            0.0,  # Vra0 - Unstressed right atrial volume
            5.0,  # K_pas_la_1 - Atrial passive curve slope, left atrium
            0.0065,  # K_pas_la_2 - Atrial passive curve exponent factor, left atrium
            0.2,  # Emax_la - Atrial active curve slope, left atrium
            0.0,  # Vla0 - Unstressed left atrial volume
            # Ventricular Model Parameters
            5.0,  # K_pas_rv_1 - Ventricular passive curve slope, right atrium
            0.003,  # K_pas_rv_2 - Ventricular passive curve exponent factor, right atrium
            0.5,  # Emax_rv - Ventricular active curve slope, right atrium
            0.0,  # Vrv0 - Unstressed right atrial volume
            2.0,  # K_pas_lv_1 - Ventricular passive curve slope, left atrium
            0.003,  # K_pas_lv_2 - Ventricular passive curve exponent factor, left atrium
            4.0,  # Emax_lv - Ventricular active curve slope, left atrium
            20.0,  # Vlv0 - Unstressed left atrial volume
            # Atrial and Ventricular Inductances and Resistances
            0.1,  # L_ra_rv - Inductance of right atrium
            10.0,  # R_ra_rv - Resistance of right atrium
            0.1,  # L_rv_pa - Inductance of right ventricle
            15.0,  # R_rv_pa - Resistance of right ventricle
            0.1,  # L_la_lv - Inductance of left atrium
            8.0,  # R_la_lv - Resistance of left atrium
            0.1,  # L_lv_ao - Inductance of left ventricle
            25.0,  # R_lv_ao - Resistance of left ventricle
            # Aortic Arch
            1000.0e-6,  # C_ao - Aortic capacitance
            # Pulmonary Resistance and Capacitance
            4000.0e-6,  # C_pa - Pulmonary capacitance
            130.0,  # R_pa - Pulmonary resistance
            # Systemic Resistance and Capacitance
            400.0e-6,  # C_sys - Systemic Capacitance
            400.0,  # R_sys_a - Systemic Resistance - Arteries
            1200.0])  # R_sys_v - Systemic Resistance - Veins

    def solve(self, tau):
        self.tau = tau
        self.t_disc = np.append(np.arange(0, self.cycleTime * self.totalCycles, tau), self.cycleTime * self.totalCycles)
        # self.t_disc = np.arange(0, self.cycleTime * self.totalCycles, tau)
        u0 = self.defIC
        g0 = np.zeros(11)
        param = self.defParam
        self.u, self.g = explicit_euler(DAE, u0, g0, param, self.t_disc)
        return self.t_disc, self.u, self.g

    def post_process(self):
        """
        :return: Array with clinical comparison data
        """
        start =int(0.5 * len(self.u))
        co_parameter = 120 / (self.cycleTime * self.totalCycles * 1000) #scale time snippet to L/min
        # HeartRate = self.defParam[0]
        RAP = np.mean(self.g[start:, 0]) * baryeTommHg # Right Atrial Pressure
        sPAP = np.max(self.u[start:, 5]) * baryeTommHg # Systolic pulmonary artery pressure
        dPAP = np.min(self.u[start:, 5]) * baryeTommHg # Diastolic pulmonary artery pressure
        PCW = np.mean(self.g[start:, 2]) * baryeTommHg # Pulmonary capillary wedge pressure - average left atrial pessure
        SBP = np.max(self.u[start:, 10]) * baryeTommHg # Systolic blood pressure
        DBP = np.min(self.u[start:, 10]) * baryeTommHg # Diastolic blood pressure
        # SVR = self.defParam[20]# Systemic vascular resistance
        CO = np.trapz(self.u[start:, 9], self.t_disc[start:]) * co_parameter # Cardial Output
        sRV = np.max(self.g[start:, 1]) * baryeTommHg # Systolic right ventricular pressure
        RVDEP = np.min(self.g[start:, 1]) * baryeTommHg # Right ventricular end diastolic pressure
        sLV = np.max(self.g[start:, 3]) * baryeTommHg # Systolic left ventricular pressure
        LVDEP = np.min(self.g[start:, 3]) * baryeTommHg # Left ventricular end diastolic pressure
        LVEF = (np.max(self.u[start:, 3]) - np.min(self.u[start:, 3])) / np.max(self.u[start:, 3]) #Left ventricular ejection fraction
        RVEF = (np.max(self.u[start:, 2]) - np.min(self.u[start:, 2])) / np.max(self.u[start:, 2])  # Right ventricular ejection fraction
        LVSV = (np.max(self.u[start:, 3]) - np.min(self.u[start:, 3])) # Left ventricular stroke volume
        RVSV = (np.max(self.u[start:, 2]) - np.min(self.u[start:, 2])) # Right ventricular stroke volume
        LVEDV = np.max(self.u[start:, 3]) # Left ventricular end diastolic volume
        RVEDV = np.max(self.u[start:, 2]) # Right ventricular end diastolic volume

        self.post=np.round([RAP, sPAP, dPAP, PCW, SBP, DBP, CO, sRV, RVDEP, sLV, LVDEP, LVEF, RVEF, LVSV, RVSV, LVEDV, RVEDV],2)


def main():

    cycleTime = 1.0
    totalCycles = 10
    tau = 1/16000

    model = HeartModel(cycleTime, totalCycles)

    # t_disc, u, g = model.solve(model.cycleTime * tau)

    # model.post_process()

    # print(model.post)

    plots.overview_plots(model, 1/1000)

    # analysis.plot_convergence(model)


    return 0


main()
