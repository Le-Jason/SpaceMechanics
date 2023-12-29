import math as m
import numpy as np

def nutationLongitudeObliquity(i,value):
    i = i - 1
    coeff_map = {"an1": 0, "an2": 1, "an3": 2, "an4": 3, "an5": 4, "A": 5, "B": 6, "C": 7, "D": 8}
    nutationCoeff = [[  0,  0,  0,  0,  1, -171996.0, -174.2,  92025.0,    8.9 ], # 1-10 
        [  0,  0,  0,  0,  2,    2062.0,    0.2,   -895.0,    0.5 ],
        [ -2,  0,  2,  0,  1,      46.0,    0.0,    -24.0,    0.0 ],
        [  2,  0, -2,  0,  0,      11.0,    0.0,      0.0,    0.0 ],
        [ -2,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ],
        [  1, -1,  0, -1,  0,      -3.0,    0.0,      0.0,    0.0 ],
        [  0, -2,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 ],
        [  2,  0, -2,  0,  1,       1.0,    0.0,      0.0,    0.0 ],
        [  0,  0,  2, -2,  2,  -13187.0,   -1.6,   5736.0,   -3.1 ],
        [  0,  1,  0,  0,  0,    1426.0,   -3.4,     54.0,   -0.1 ],
        [  0,  1,  2, -2,  2,    -517.0,    1.2,    224.0,   -0.6 ], # 11-20 
        [  0, -1,  2, -2,  2,     217.0,   -0.5,    -95.0,    0.3 ],
        [  0,  0,  2, -2,  1,     129.0,    0.1,    -70.0,    0.0 ],
        [  2,  0,  0, -2,  0,      48.0,    0.0,      1.0,    0.0 ],
        [  0,  0,  2, -2,  0,     -22.0,    0.0,      0.0,    0.0 ],
        [  0,  2,  0,  0,  0,      17.0,   -0.1,      0.0,    0.0 ],
        [  0,  1,  0,  0,  1,     -15.0,    0.0,      9.0,    0.0 ],
        [  0,  2,  2, -2,  2,     -16.0,    0.1,      7.0,    0.0 ],
        [  0, -1,  0,  0,  1,     -12.0,    0.0,      6.0,    0.0 ],
        [ -2,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 ],
        [  0, -1,  2, -2,  1,      -5.0,    0.0,      3.0,    0.0 ], # 21-30 
        [  2,  0,  0, -2,  1,       4.0,    0.0,     -2.0,    0.0 ],
        [  0,  1,  2, -2,  1,       4.0,    0.0,     -2.0,    0.0 ],
        [  1,  0,  0, -1,  0,      -4.0,    0.0,      0.0,    0.0 ],
        [  2,  1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 ],
        [  0,  0, -2,  2,  1,       1.0,    0.0,      0.0,    0.0 ],
        [  0,  1, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [  0,  1,  0,  0,  2,       1.0,    0.0,      0.0,    0.0 ],
        [ -1,  0,  0,  1,  1,       1.0,    0.0,      0.0,    0.0 ],
        [  0,  1,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [  0,  0,  2,  0,  2,   -2274.0,   -0.2,    977.0,   -0.5 ], # 31-40 
        [  1,  0,  0,  0,  0,     712.0,    0.1,     -7.0,    0.0 ],
        [  0,  0,  2,  0,  1,    -386.0,   -0.4,    200.0,    0.0 ],
        [  1,  0,  2,  0,  2,    -301.0,    0.0,    129.0,   -0.1 ],
        [  1,  0,  0, -2,  0,    -158.0,    0.0,     -1.0,    0.0 ],
        [ -1,  0,  2,  0,  2,     123.0,    0.0,    -53.0,    0.0 ],
        [  0,  0,  0,  2,  0,      63.0,    0.0,     -2.0,    0.0 ],
        [  1,  0,  0,  0,  1,      63.0,    0.1,    -33.0,    0.0 ],
        [ -1,  0,  0,  0,  1,     -58.0,   -0.1,     32.0,    0.0 ],
        [ -1,  0,  2,  2,  2,     -59.0,    0.0,     26.0,    0.0 ],
        [  1,  0,  2,  0,  1,     -51.0,    0.0,     27.0,    0.0 ], # 41-50 
        [  0,  0,  2,  2,  2,     -38.0,    0.0,     16.0,    0.0 ],
        [  2,  0,  0,  0,  0,      29.0,    0.0,     -1.0,    0.0 ],
        [  1,  0,  2, -2,  2,      29.0,    0.0,    -12.0,    0.0 ],
        [  2,  0,  2,  0,  2,     -31.0,    0.0,     13.0,    0.0 ],
        [  0,  0,  2,  0,  0,      26.0,    0.0,     -1.0,    0.0 ],
        [ -1,  0,  2,  0,  1,      21.0,    0.0,    -10.0,    0.0 ],
        [ -1,  0,  0,  2,  1,      16.0,    0.0,     -8.0,    0.0 ],
        [  1,  0,  0, -2,  1,     -13.0,    0.0,      7.0,    0.0 ],
        [ -1,  0,  2,  2,  1,     -10.0,    0.0,      5.0,    0.0 ],
        [  1,  1,  0, -2,  0,      -7.0,    0.0,      0.0,    0.0 ], # 51-60 
        [  0,  1,  2,  0,  2,       7.0,    0.0,     -3.0,    0.0 ],
        [  0, -1,  2,  0,  2,      -7.0,    0.0,      3.0,    0.0 ],
        [  1,  0,  2,  2,  2,      -8.0,    0.0,      3.0,    0.0 ],
        [  1,  0,  0,  2,  0,       6.0,    0.0,      0.0,    0.0 ],
        [  2,  0,  2, -2,  2,       6.0,    0.0,     -3.0,    0.0 ],
        [  0,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 ],
        [  0,  0,  2,  2,  1,      -7.0,    0.0,      3.0,    0.0 ],
        [  1,  0,  2, -2,  1,       6.0,    0.0,     -3.0,    0.0 ],
        [  0,  0,  0, -2,  1,      -5.0,    0.0,      3.0,    0.0 ],
        [  1, -1,  0,  0,  0,       5.0,    0.0,      0.0,    0.0 ], # 61-70 
        [  2,  0,  2,  0,  1,      -5.0,    0.0,      3.0,    0.0 ],
        [  0,  1,  0, -2,  0,      -4.0,    0.0,      0.0,    0.0 ],
        [  1,  0, -2,  0,  0,       4.0,    0.0,      0.0,    0.0 ],
        [  0,  0,  0,  1,  0,      -4.0,    0.0,      0.0,    0.0 ],
        [  1,  1,  0,  0,  0,      -3.0,    0.0,      0.0,    0.0 ],
        [  1,  0,  2,  0,  0,       3.0,    0.0,      0.0,    0.0 ],
        [  1, -1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ],
        [ -1, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 ],
        [ -2,  0,  0,  0,  1,      -2.0,    0.0,      1.0,    0.0 ],
        [  3,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ], # 71-80 
        [  0, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 ],
        [  1,  1,  2,  0,  2,       2.0,    0.0,     -1.0,    0.0 ],
        [ -1,  0,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 ],
        [  2,  0,  0,  0,  1,       2.0,    0.0,     -1.0,    0.0 ],
        [  1,  0,  0,  0,  2,      -2.0,    0.0,      1.0,    0.0 ],
        [  3,  0,  0,  0,  0,       2.0,    0.0,      0.0,    0.0 ],
        [  0,  0,  2,  1,  2,       2.0,    0.0,     -1.0,    0.0 ],
        [ -1,  0,  0,  0,  2,       1.0,    0.0,     -1.0,    0.0 ],
        [  1,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [ -2,  0,  2,  2,  2,       1.0,    0.0,     -1.0,    0.0 ], # 81-90 
        [ -1,  0,  2,  4,  2,      -2.0,    0.0,      1.0,    0.0 ],
        [  2,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [  1,  1,  2, -2,  2,       1.0,    0.0,     -1.0,    0.0 ],
        [  1,  0,  2,  2,  1,      -1.0,    0.0,      1.0,    0.0 ],
        [ -2,  0,  2,  4,  2,      -1.0,    0.0,      1.0,    0.0 ],
        [ -1,  0,  4,  0,  2,       1.0,    0.0,      0.0,    0.0 ],
        [  1, -1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 ],
        [  2,  0,  2, -2,  1,       1.0,    0.0,     -1.0,    0.0 ],
        [  2,  0,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0 ],
        [  1,  0,  0,  2,  1,      -1.0,    0.0,      0.0,    0.0 ], # 91-100
        [  0,  0,  4, -2,  2,       1.0,    0.0,      0.0,    0.0 ],
        [  3,  0,  2, -2,  2,       1.0,    0.0,      0.0,    0.0 ],
        [  1,  0,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [  0,  1,  2,  0,  1,       1.0,    0.0,      0.0,    0.0 ],
        [ -1, -1,  0,  2,  1,       1.0,    0.0,      0.0,    0.0 ],
        [  0,  0, -2,  0,  1,      -1.0,    0.0,      0.0,    0.0 ],
        [  0,  0,  2, -1,  2,      -1.0,    0.0,      0.0,    0.0 ],
        [  0,  1,  0,  2,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [  1,  0, -2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [  0, -1,  2,  0,  1,      -1.0,    0.0,      0.0,    0.0 ], # 101-106 
        [  1,  1,  0, -2,  1,      -1.0,    0.0,      0.0,    0.0 ],
        [  1,  0, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 ],
        [  2,  0,  0,  2,  0,       1.0,    0.0,      0.0,    0.0 ],
        [  0,  0,  2,  4,  2,      -1.0,    0.0,      0.0,    0.0 ],
        [  0,  1,  0,  1,  0,       1.0,    0.0,      0.0,    0.0 ]
    ]
    return nutationCoeff[i][coeff_map[value]]

def ROT1(alpha):
    A = np.array([1.0, 0.0, 0.0,
                    0.0, m.cos(alpha), m.sin(alpha),
                    0.0, -m.sin(alpha), m.cos(alpha)])
    A = A.reshape((3,3))
    return A

def ROT2(alpha):
    A = np.array([m.cos(alpha), 0.0, -m.sin(alpha),
                    0.0, 1.0, 0.0,
                    m.sin(alpha), 0.0, m.cos(alpha)])
    A = A.reshape((3,3))
    return A

def ROT3(alpha):
    A = np.array([m.cos(alpha), m.sin(alpha), 0.0,
                    -m.sin(alpha), m.cos(alpha), 0.0,
                    0.0, 0.0, 1.0])
    A = A.reshape((3,3))
    return A

def format_time(seconds):
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    return f"{hours:02f}:{minutes:02f}:{seconds:02f}"

def DMS2Rad(deg,arcmin,arcsec):
    alpha = (deg + (arcmin/60) + (arcsec/3600)) * (np.pi/180)
    return alpha

def Rad2DMS(alpha):
    temp = alpha * (180/np.pi)
    deg = int(temp)
    arcmin = int((temp - deg)*60)
    arcsec = (temp - deg - (arcmin/60))*3600
    return deg,arcmin,arcsec

def HMS2Rad(h,min,s):
    t = 15*(h + (min/60) + (s/3600))*(np.pi/180)
    return t

def Rad2HMS(t):
    temp = t*(180/(15*np.pi))
    h = int(temp)
    min = int((temp-h)*60)
    s = (temp - h - (min/60))*3600
    return h,min,s

def julianDate(yr,mo,d,h,min,s):
    JD = (367*yr) - int((7*(yr + int((mo+9)/12)))/4) + int((275*mo/9)) + d + 1721013.5 + (((((s/60)+min)/60)+h)/24)
    return JD

def convTime(yr,mo,day,UTC,delta_UT1,delta_AT):
    UT1 = UTC + delta_UT1
    TAI = UTC + delta_AT
    GPS = UTC + delta_AT - 19
    TT = TAI + 32.184
    TT_hr = TT//3600
    TT_min = (TT % 3600) // 60
    TT_sec = TT % 60
    JD_TT = julianDate(yr,mo,day,TT_hr,TT_min,TT_sec)
    T_TT = (JD_TT - 2451545.0)/36525.0
    TAI_hr = TAI//3600
    TAI_min = (TAI % 3600) // 60
    TAI_sec = TAI % 60
    JD_TAI = julianDate(yr,mo,day,TAI_hr,TAI_min,TAI_sec)
    UT1_hr = UT1//3600
    UT1_min = (UT1 % 3600) // 60
    UT1_sec = UT1 % 60
    JD_UT1 = julianDate(yr,mo,day,UT1_hr,UT1_min,UT1_sec)
    T_UT1 = (JD_UT1 - 2451545.0)/36525.0
    # TDB is slightly off
    # This TDB needs more terms for more accuracy
    TDB = (TT + (0.001657)*m.sin(628.3076*T_TT + 6.2401)
            + 0.000022*m.sin(575.3385*T_TT + 4.2970) + 0.000014*m.sin(1256.6152*T_TT + 6.1969)
            + 0.000005*m.sin(606.9777*T_TT + 4.0212) + 0.000005*m.sin(52.9691*T_TT + 0.4444)
            + 0.000002*m.sin(21.3299*T_TT + 5.5431) + 0.000010*T_TT*m.sin(628.3076*T_TT + 4.2490))
    TDB_hr = TAI//3600
    TDB_min = (TAI % 3600) // 60
    TDB_sec = TAI % 60
    JD_TDB = julianDate(yr,mo,day,TDB_hr,TDB_min,TDB_sec)
    T_TDB = (JD_TDB - 2451545.0)/36525.0
    return UT1,TAI,TT,TDB,T_UT1,T_TT,T_TDB


def FK5(r_GCRF,v_GCRF,yr,mo,day,UTC,delta_UT1,delta_AT,xp,yp,delta_delta_psi_1980=0,delta_delta_epsilon_1980=0,LOD=0,itrf2gcrf=True):
    #Using 1980 and 1982 values

    w = 7.29211514670698e-05 * (1.0  - LOD/86400.0)
    w_earth = np.array([0,0,w])
    w_earth = w_earth.reshape((3,1))

    # These values are in seconds units
    UT1,TAI,TT,TDB,T_UT1,T_TT,T_TDB = convTime(yr,mo,day,UTC,delta_UT1,delta_AT)

    # These values are in arcsec units
    zeta = (2306.2181*T_TT) + (0.30188*(T_TT**2)) + (0.017998*(T_TT**3))
    theta = (2004.3109*T_TT) - (0.42665*(T_TT**2)) - (0.041833*(T_TT**3))
    z = (2306.2181*T_TT) + (1.09468*(T_TT**2)) + (0.018203*(T_TT**3))
    zeta = zeta * (1/3600) * (np.pi/180)
    theta = theta * (1/3600) * (np.pi/180)
    z = z * (1/3600) * (np.pi/180)

    # This value is in arcsec
    epsilon_line_1980 = 84381.448 - (46.8150*T_TT) - (0.00059*(T_TT**2)) + (0.001813*(T_TT**3))

    # These values are in degs
    M_moon = m.fmod(((((  0.064 ) * T_TT + 31.310 ) * T_TT + 1717915922.6330 ) * T_TT ) / 3600.0 + 134.96298139,360) * (np.pi/180)
    M_sun = m.fmod(((((- 0.012 ) * T_TT -  0.577 ) * T_TT +  129596581.2240 ) * T_TT ) / 3600.0 + 357.52772333,360) * (np.pi/180)
    u_M_moon = m.fmod(((((  0.011 ) * T_TT - 13.257 ) * T_TT + 1739527263.1370 ) * T_TT ) / 3600.0 +  93.27191028,360) * (np.pi/180)
    D_sun = m.fmod(((((  0.019 ) * T_TT -  6.891 ) * T_TT + 1602961601.3280 ) * T_TT ) / 3600.0 + 297.85036306,360) * (np.pi/180)
    omega_moon = m.fmod(((((  0.008 ) * T_TT +  7.455 ) * T_TT - 6962890.5390 ) * T_TT ) / 3600.0 + 125.04452222,360) * (np.pi/180)

    #Finding the nutation in obliquity
    # NOTE: delta_delta_psi_1980 & delta_delta_epsilon_1980 are in deg
    delta_epsilon_1980 = 0
    delta_psi_1980 = 0
    convert = 0.0001 / 3600.0 # 0.0001" to deg
    for i in range(1,107):
        api = ((nutationLongitudeObliquity(i,"an1")*M_moon) +
                (nutationLongitudeObliquity(i,"an2")*M_sun) +
                (nutationLongitudeObliquity(i,"an3")*u_M_moon) + 
                (nutationLongitudeObliquity(i,"an4")*D_sun) +
                (nutationLongitudeObliquity(i,"an5")*omega_moon))
        delta_epsilon_1980 = delta_epsilon_1980 + ((nutationLongitudeObliquity(i,"C")*convert) + ((nutationLongitudeObliquity(i,"D")*convert)*T_TT)) * m.cos(api)
        delta_psi_1980 = delta_psi_1980 + ((nutationLongitudeObliquity(i,"A")*convert) + ((nutationLongitudeObliquity(i,"B")*convert)*T_TT)) * m.sin(api)
    delta_epsilon_1980 = m.fmod(delta_epsilon_1980 + delta_delta_epsilon_1980,360.0)
    delta_psi_1980 = m.fmod(delta_psi_1980 + delta_delta_psi_1980,360.0)
    epsilon_1980 = (delta_epsilon_1980) + (epsilon_line_1980/3600.0)
    epsilon_line_1980 = (epsilon_line_1980/3600.0) * (np.pi/180)
    epsilon_1980 = epsilon_1980 * (np.pi/180)
    delta_psi_1980 = delta_psi_1980 * (np.pi/180)

    # GMST was in seconds
    theta_GMST = 67310.54841 + ((876600.0 * 3600 + 8640184.812866)*T_UT1) +( 0.093104*(T_UT1**2)) - ((6.2e-6)*(T_UT1**3))
    theta_GMST = m.fmod(theta_GMST * (np.pi/180) / 240.0, 2*np.pi)
    if theta_GMST < 0.0:
        theta_GMST += (2*np.pi)
    theta_GMST = theta_GMST * (180/np.pi)
    if (True): # UT1 > 2450449.5
        Eq_equniox_1982 = (delta_psi_1980*(180/np.pi)*m.cos(epsilon_1980)) + (0.00264*(1/3600)*m.sin(omega_moon)) + (0.000063*(1/3600)*m.sin(2*omega_moon))
    else:
        Eq_equniox_1982 = (delta_psi_1980*(180/np.pi)*m.cos(epsilon_1980))
    theta_GAST_1982 = theta_GMST + Eq_equniox_1982
    theta_GAST_1982 = m.fmod(theta_GAST_1982,360) * (np.pi/180)
    
    # 3-2-3 Rotation ROT3(zeta)*ROT2(-theta)*ROT3(z)
    P_Matrix = np.matmul(ROT2(-theta), ROT3(z))
    P_Matrix = np.matmul(ROT3(zeta), P_Matrix)

    # 1-3-1 Rotation ROT1(-epsilon_line_1980)*ROT3(delta_psi_1980)*ROT1(epsilon_1980)
    N_Matrix = np.matmul(ROT3(delta_psi_1980), ROT1(epsilon_1980))
    N_Matrix = np.matmul(ROT1(-epsilon_line_1980), N_Matrix)
    
    R_Matrix = ROT3(-theta_GAST_1982)

    # 1-2 Rotation ROT1(yp)*ROT2(xp)
    W_Matrix = np.matmul(ROT1(yp*(np.pi/180)), ROT2(xp*(np.pi/180)))

    if itrf2gcrf:
        # r_ITRF = P * N * R * W * r_GCRF
        T_GCRF_ITRF = np.matmul(P_Matrix, N_Matrix)
        T_GCRF_ITRF = np.matmul(T_GCRF_ITRF, R_Matrix)
        T_GCRF_ITRF = np.matmul(T_GCRF_ITRF, W_Matrix)
        r_ITRF = np.matmul(T_GCRF_ITRF, r_GCRF)

        # r_PEF = ROT1(yp) * ROT2(xp) * r_ITRF
        T_PEF_ITRF = np.matmul(ROT1(yp*(np.pi/180)), ROT2(xp*(np.pi/180)))
        r_PEF = np.matmul(T_PEF_ITRF, r_GCRF)

        # v_ITRF = P * N * R * (W * v_GCRF - w_earth x r_PEF)
        T_PEF_ITRF_vel = np.matmul(P_Matrix, N_Matrix)
        T_PEF_ITRF_vel = np.matmul(T_PEF_ITRF_vel, R_Matrix)
        v_ITRF = np.matmul(W_Matrix, v_GCRF)
        cross_pro = np.cross(w_earth.flatten(),r_PEF.flatten())
        cross_pro = cross_pro.reshape((3,1))
        v_ITRF = v_ITRF + cross_pro
        v_ITRF = np.matmul(T_PEF_ITRF_vel,v_ITRF)
    else:
        # r_ITRF = W_T * R_T * N_T * P_T * r_GCRF
        T_ITRF_GCRF = np.matmul(W_Matrix.T, R_Matrix.T)
        T_ITRF_GCRF = np.matmul(T_ITRF_GCRF, N_Matrix.T)
        T_ITRF_GCRF = np.matmul(T_ITRF_GCRF, P_Matrix.T)
        r_ITRF = np.matmul(T_ITRF_GCRF, r_GCRF)

        # r_PEF = ROT1(yp) * ROT2(xp) * r_ITRF
        T_PEF_ITRF = np.matmul(ROT1(yp*(np.pi/180)), ROT2(xp*(np.pi/180)))
        r_PEF = np.matmul(T_PEF_ITRF, r_ITRF)

        # v_ITRF = W_T * (R_T * N_T * P_T * v_GCRF - w_earth x r_PEF)
        T_PEF_ITRF_vel = np.matmul(R_Matrix.T, N_Matrix.T)
        T_PEF_ITRF_vel = np.matmul(T_PEF_ITRF_vel, P_Matrix.T)
        v_ITRF = np.matmul(T_PEF_ITRF_vel, v_GCRF)
        cross_pro = np.cross(w_earth.flatten(),r_PEF.flatten())
        cross_pro = cross_pro.reshape((3,1))
        v_ITRF = v_ITRF - cross_pro
        v_ITRF = np.matmul(W_Matrix.T,v_ITRF)

    return r_ITRF, v_ITRF


def siteTrack(phi_gd,lambda_,h_ellp,rho,B,el,rho_dot,B_dot,el_dot,yr,mo,day,UTC,delta_UT1,delta_AT,xp,yp):
    R_earth = 6378.1363
    e_earth = 0.081819221456
    C_earth = R_earth / m.sqrt(1 - ((e_earth**2)*(m.sin(phi_gd)**2)))
    S_earth = C_earth*(1-e_earth**2)
    r_delta = (C_earth + h_ellp)*m.cos(phi_gd)
    r_k = (S_earth + h_ellp)*m.sin(phi_gd)
    r_site_ECEF = [r_delta*m.cos(lambda_),r_delta*m.sin(lambda_),r_k]
    r_site_ECEF = np.array(r_site_ECEF)
    r_site_ECEF = r_site_ECEF.reshape((3,1))
    rho_SEZ = [-rho*m.cos(el)*m.cos(B),rho*m.cos(el)*m.sin(B),rho*m.sin(el)]
    rho_SEZ = np.array(rho_SEZ)
    rho_SEZ = rho_SEZ.reshape((3,1))
    rho_dot_SEZ = [(-rho_dot*m.cos(el)*m.cos(B)) + (rho*m.sin(el)*m.cos(B)*el_dot) + (rho*m.cos(el)*m.sin(B)*B_dot),
                    (rho_dot*m.cos(el)*m.sin(B)) - (rho*m.sin(el)*m.sin(B)*el_dot) + (rho*m.cos(el)*m.cos(B)*B_dot),
                    (rho_dot*m.sin(el)) + (rho*m.cos(el)*el_dot)]
    rho_dot_SEZ = np.array(rho_dot_SEZ)
    rho_dot_SEZ = rho_dot_SEZ.reshape((3,1))
    T_ECEF_SEZ = np.matmul(ROT3(-lambda_), ROT2(-(np.deg2rad(90)-phi_gd)))
    rho_ECEF = np.matmul(T_ECEF_SEZ,rho_SEZ)
    rho_dot_ECEF = np.matmul(T_ECEF_SEZ,rho_dot_SEZ)
    r_ECEF = rho_ECEF + r_site_ECEF
    v_ECEF = rho_dot_ECEF
    r_ECI,v_ECI = FK5(r_ECEF,v_ECEF,yr,mo,day,UTC,delta_UT1,delta_AT,xp,yp)
    return r_ECI,v_ECI

def anglesOnlyGauss(L1,L2,L3,JD1,JD2,JD3,r_site1,r_site2,r_site3):
    return

deg2rad = np.pi/180
r = 6378.0

r_ITRF = np.array([-1033.4793830, 7901.2952754, 6380.3565958])
v_ITRF = np.array([-3.225636520, -2.872451450, 5.531924446])
r_ITRF = r_ITRF.reshape((3,1))
v_ITRF = v_ITRF.reshape((3,1))

#print(FK5(r_ITRF,v_ITRF,2004,4,6,(7*3600) + (51*60) + (28.386009),-0.4399619,32.0,-0.140682*(1/3600),0.333309*(1/3600),delta_delta_psi_1980=-0.052195*(1/3600),delta_delta_epsilon_1980=-0.003875*(1/3600),LOD=0.0015563))
#print(siteTrack(39.007*deg2rad,-104.883*deg2rad,2187/1000,604.68,205.6*deg2rad,30.7*deg2rad,2.08,0.15*deg2rad,0.17*deg2rad,1995,5,20,3*3600 + 17*60 +2, 0, 29, 0.0, 0.0))