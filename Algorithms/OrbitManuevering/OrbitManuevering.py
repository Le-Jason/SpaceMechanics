import math as m
import numpy as np

def OneTangentBurn(r_i,r_f,trueanomaly_trans_b):
    k = 0
    E_0 = 0
    mu = 398600.4418
    R_neg1 = r_i/r_f
    if trueanomaly_trans_b <= np.deg2rad(180):
        e_trans = (R_neg1 - 1) / (m.cos(trueanomaly_trans_b) - R_neg1)
        a_trans = r_i/(1 - e_trans)
    else:
        e_trans = (R_neg1 - 1) / (m.cos(trueanomaly_trans_b) + R_neg1)
        a_trans = r_i/(1 + e_trans)
    v_i = m.sqrt(mu/r_i)
    v_f = m.sqrt(mu/r_f)
    v_trans_a = m.sqrt((2*mu/r_i) - (mu/a_trans))
    v_trans_b = m.sqrt((2*mu/r_f) - (mu/a_trans))
    delta_v_a = v_trans_a - v_i
    phi_trans_b = m.atan2((e_trans*m.sin(trueanomaly_trans_b)),(1 + (e_trans*m.cos(trueanomaly_trans_b))))
    print(trueanomaly_trans_b*180/np.pi)
    print(phi_trans_b*180/np.pi)
    delta_v_b = m.sqrt((v_trans_b**2) + (v_f**2) - (2*v_trans_b*v_f*m.cos(phi_trans_b)))
    delta_v_otb = abs(delta_v_a) + abs(delta_v_b)
    E = m.acos((e_trans + m.cos(trueanomaly_trans_b))/(1 + (e_trans*m.cos(trueanomaly_trans_b))))
    period_trans = m.sqrt((a_trans**3)/mu) * ((2*k*np.pi) + (E - e_trans*m.sin(E)) - (E_0 - e_trans*m.sin(E_0)))
    return e_trans,a_trans,period_trans,delta_v_a,delta_v_b

def inclinationOnlyCircular(delta_i,phi_fpa,v_i):
    #circular orbits only
    delta_v_ionly = 2*v_i*m.cos(phi_fpa)*m.sin(delta_i/2)
    return delta_v_ionly

def inclinationOnlyElliptical(delta_i,e,p,w,trueAnomaly):
    #Elliptical orbits only
    mu = 398600.4418
    a = p / (1 - e**2)
    r = p / (1 + e*m.cos(trueAnomaly))
    v = m.sqrt((2*mu/r) - (mu/a))
    phi_fpa = m.atan2(e*m.sin(trueAnomaly),1 + e*m.cos(trueAnomaly))
    delta_v_ionly_1 = 2*v*m.cos(phi_fpa)*m.sin(delta_i/2)

    #check other node point
    if np.rad2deg(trueAnomaly) < 180:
        trueAnomaly = np.deg2rad(np.rad2deg(trueAnomaly) + 180)
    else:
        trueAnomaly = np.deg2rad(np.rad2deg(trueAnomaly) - 180)
    r = p / (1 + e*m.cos(trueAnomaly))
    v = m.sqrt((2*mu/r) - (mu/a))
    phi_fpa = m.atan2(e*m.sin(trueAnomaly),1 + e*m.cos(trueAnomaly))
    delta_v_ionly_2 = 2*v*m.cos(phi_fpa)*m.sin(delta_i/2)

    return delta_v_ionly_1,delta_v_ionly_2

def ascendingNodeCircular(delta_RAAN,i_i,v_i):
    #Circular orbits only
    angle = m.acos((m.cos(i_i)**2) + (m.sin(i_i)**2)*m.cos(delta_RAAN))
    delta_v_RAAN = 2*v_i*m.sin(angle/2)
    return delta_v_RAAN

def inclinationAscendingNodeCircular(i_i,i_f,delta_RAAN,v_i):
    #Circular Orbits Only
    angle = m.acos((m.cos(i_i)*m.cos(i_f))+(m.sin(i_i)*m.sin(i_f)*m.cos(delta_RAAN)))
    delta_v = 2*v_i*m.sin(angle/2)
    return delta_v

def minCombinedPlaneChangeCircular(i_i,r_i,v_i,v_trans_a,v_trans_b,r_f,v_f,delta_i):
    #Circular Orbits Only
    iterate = 0
    if iterate == 1:
        s = 0
        part1 = m.sin(s*delta_i)
        delta_v_a = m.sqrt((v_i**2) + (v_trans_a**2) - (2*v_i*v_trans_a*m.cos(s*delta_i)))
        delta_v_b = m.sqrt((v_f**2) + (v_trans_b**2) - (2*v_f*v_trans_b*m.cos((1-s)*delta_i)))
        part2 = (delta_v_a*v_f*v_trans_b*m.sin((1-s)*delta_i))/(delta_v_b*v_i*v_trans_a)
        tol = 0.0001
        while ((part1 + tol) >= part2) and ((part1 - tol) <= part2):
            s += 0.00001
            part1 = m.sin(s*delta_i)
            delta_v_a = m.sqrt((v_i**2) + (v_trans_a**2) - (2*v_i*v_trans_a*m.cos(s*delta_i)))
            delta_v_b = m.sqrt((v_f**2) + (v_trans_b**2) - (2*v_f*v_trans_b*m.cos((1-s)*delta_i)))
            part2 = (delta_v_a*v_f*v_trans_b*m.sin((1-s)*delta_i))/(delta_v_b*v_i*v_trans_a)
    else:
        R = r_f/r_i
        s = (1/delta_i)*m.atan2(m.sin(delta_i),(R**(3/2) + m.cos(delta_i)))
        delta_v_a = m.sqrt((v_i**2) + (v_trans_a**2) - (2*v_i*v_trans_a*m.cos(s*delta_i)))
        delta_v_b = m.sqrt((v_f**2) + (v_trans_b**2) - (2*v_f*v_trans_b*m.cos((1-s)*delta_i)))
    delta_i_i = s*delta_i
    delta_i_f = (1-s)*delta_i
    delta_v_i = m.sqrt((v_i**2) + (v_trans_a**2) - (2*v_i*v_trans_a*m.cos(delta_i_i)))
    delta_v_f = m.sqrt((v_f**2) + (v_trans_b**2) - (2*v_f*v_trans_b*m.cos(delta_i_f)))

    return delta_i_i,delta_i_f,delta_v_a,delta_v_b

def fixDeltaVManeuvers(v_i,v_f,v_trans_a,delta_v):
    gramma = m.acos(-1*((v_i**2) + (delta_v**2) - (v_f**2))/(2*v_i*delta_v))
    delta_i_1 = ((v_i**2) + (v_trans_a**2) - (delta_v**2))/(2*v_i*v_trans_a)
    return gramma,delta_i_1

def circularCoplanerPhasingSame(a_tgt,phase_ang,k_tgt,k_int):
    #Same Orbits
    mu = 398600.4418
    w_tgt = m.sqrt(mu/a_tgt**3)
    period_phase = (2*np.pi*k_tgt + phase_ang)/w_tgt
    a_phase = (mu*(period_phase/(2*np.pi*k_int))**2)**(1/3)
    delta_v = 2*abs(m.sqrt((2*mu/a_tgt) - (mu/a_phase)) - (m.sqrt(mu/a_tgt)))
    return period_phase,delta_v,a_phase

def circularCoplanerPhasingDifferent(phase_ang_i,a_int,a_tgt,k):
    #Different Orbits
    mu = 398600.4418
    w_tgt = m.sqrt(mu/a_tgt**3)
    w_int = m.sqrt(mu/a_int**3)
    a_trans = (a_int + a_tgt)/2
    period_trans = np.pi*m.sqrt((a_trans**3)/mu)
    alpha_L = w_tgt*period_trans
    phase_ang = alpha_L - np.pi
    period_wait = (phase_ang - phase_ang_i + (2*np.pi*k))/(w_int - w_tgt)
    delta_v = abs(m.sqrt((2*mu/a_int) - (mu/a_trans)) - m.sqrt(mu/a_int)) + abs(m.sqrt((2*mu/a_tgt) - (mu/a_trans)) - m.sqrt(mu/a_tgt))
    return period_trans,delta_v,a_trans

def nonCoplanarPhasing(phase_ang_i,a_int,a_tgt,k_tgt,u_int,RAAN_int,lambda_true_0,delta_i):
    #Circular Orbits Only
    mu = 398600.4418
    k_int = 1 #assume
    w_tgt = m.sqrt(mu/(a_tgt**3))
    w_int = m.sqrt(mu/(a_int**3))
    a_trans = (a_int + a_tgt)/2
    period_trans = np.pi*m.sqrt((a_trans**3)/mu)
    alpha_L = w_tgt*period_trans
    #Find delta_phase_ang_int to reach a node (180 or 360 - u_int)
    delta_phase_ang_int = np.deg2rad(180) - u_int
    delta_t_node = delta_phase_ang_int/w_int
    lambda_true_tgt_1 = lambda_true_0 + (w_tgt*delta_t_node)
    #Find lambda_true for interceptor at t_1 lambda_true_int_1 = RAAN + pi
    lambda_true_int_1 = RAAN_int + np.pi
    phase_ang_new = lambda_true_int_1 - lambda_true_tgt_1
    alpha_new = np.pi + phase_ang_new
    period_phase = (alpha_new - alpha_L + (2*np.pi*k_tgt))/w_tgt
    a_phase = (mu*(period_phase/(k_int*2*np.pi))**2)**(1/3)
    delta_v_phase = abs(m.sqrt((2*mu/a_int) - (mu/a_phase)) - m.sqrt(mu/a_int))
    delta_v_trans_1 = abs(m.sqrt((2*mu/a_int) - (mu/a_trans)) - m.sqrt((2*mu/a_int) - (mu/a_phase)))
    delta_v_trans_2 = m.sqrt(((2*mu/a_tgt) - (mu/a_trans)) + (mu/a_tgt) - (2*m.sqrt((2*mu/a_tgt) - (mu/a_trans))*m.sqrt(mu/a_tgt)*m.cos(delta_i)))
    period_total = 2*np.pi*m.sqrt((a_phase**3)/mu) + period_trans + delta_t_node
    delta_v_tot = delta_v_phase + delta_v_trans_1 + delta_v_trans_2

    return period_trans,period_phase,delta_v_phase,delta_v_trans_1,delta_v_trans_2,a_phase

# def lowThrustTransfer(a_i,a_f,i_i,i_f,Isp,F,m_i):
#     g = 9.80665
#     mu = 398600.4418
#     R = a_f/a_i
#     m_dot = -F/(Isp*g)
#     m_dot_spec = m_dot/m_i
#     a_thrust_i = F/m_i
#     delta_i = i_f - i_i
#     theta_p = 0.0
#     delta_t = 0.0
#     r0 = a_i 
#     TU = m.sqrt((r0**3)/mu)
#     c_v = 1/(4*(lambdai**2)*(a_curr**2) + 1)
#     phi_c = m.atan2(m.cos(theta_p),m.sqrt((1/c_v) - 1))
#     a_thrust = a_thrust_i/(1 + (m_dot_spec*delta_t))
#     a_NTW = a_thrust*[0,m.cos(phi_c),m.sin(phi_c)]
#     a_NTW = np.array(a_NTW)
#     a_NTW = a_NTW.reshape((3,1))
#     R_ECI_NTW = 
#     a_ECI = np.matmul(R_ECI_NTW,a_NTW)

def hillEquations(x_0,y_0,z_0,x_dot_0,y_dot_0,z_dot_0,w_tgt,delta_t):
    #nearly circular orbit
    x = ((x_dot_0/w_tgt)*m.sin(w_tgt*delta_t)) - (((3*x_0) + (2*y_dot_0/w_tgt))*m.cos(w_tgt*delta_t)) + ((4*x_0) + (2*y_dot_0/w_tgt))
    y = (((6*x_0) + (4*y_dot_0/w_tgt))*m.sin(w_tgt*delta_t)) + (2*x_dot_0*m.cos(w_tgt*delta_t)/w_tgt) - (((6*w_tgt*x_0) + (3*y_dot_0))*delta_t) + (y_0 - (2*x_dot_0/w_tgt))
    z = (z_0*m.cos(w_tgt*delta_t)) + (z_dot_0*m.sin(w_tgt*delta_t)/w_tgt)
    x_dot = (x_dot_0*m.cos(w_tgt*delta_t)) + ((3*w_tgt*x_0) + (2*y_dot_0))*m.sin(w_tgt*delta_t)
    y_dot = (((6*w_tgt*x_0) + (4*y_dot_0))*m.cos(w_tgt*delta_t)) - (2*x_dot_0*m.sin(w_tgt*delta_t)) - ((6*w_tgt*x_0) + (3*y_dot_0))
    z_dot = (-z_0*w_tgt*m.sin(w_tgt*delta_t)) + (z_dot_0*m.cos(w_tgt*delta_t))
    return x,y,z,x_dot,y_dot,z_dot

# def hillEQCM2ECI()
    
# def ECI2HillEQCM()

deg2rad = np.pi/180
r = 6378.0
print(hillEquations(0,0,0,-0.1,-0.04,-0.02,0.0010854,20*60))