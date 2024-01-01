import math as m
import numpy as np
from scipy.optimize import fsolve
import sys
import os
current_dir = os.getcwd()
full_path = os.path.join(current_dir, "Algorithms", "KeplerProblems")
sys.path.append(full_path)  # Replace '/path/to/directory' with the actual path
from KeplerProblems import RV2COE
from sympy import symbols, Eq, solve

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
    return UT1,TAI,TT,TDB,T_UT1,T_TT,T_TDB,JD_TT,JD_TAI,JD_UT1,JD_TDB


def FK5(r_GCRF,v_GCRF,yr,mo,day,UTC,delta_UT1,delta_AT,xp,yp,delta_delta_psi_1980=0,delta_delta_epsilon_1980=0,LOD=0,itrf2gcrf=True):
    #Using 1980 and 1982 values

    w = 7.29211514670698e-05 * (1.0  - LOD/86400.0)
    w_earth = np.array([0,0,w])
    w_earth = w_earth.reshape((3,1))

    # These values are in seconds units
    UT1,TAI,TT,TDB,T_UT1,T_TT,T_TDB,JD_TT,JD_TAI,JD_UT1,JD_TDB = convTime(yr,mo,day,UTC,delta_UT1,delta_AT)

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

def Gibbs(r1,r2,r3):
    mu = 398600.4418
    Z12 = np.cross(r1.flatten(),r2.flatten())
    Z23 = np.cross(r2.flatten(),r3.flatten())
    Z31 = np.cross(r3.flatten(),r1.flatten())
    
    # Must be low
    alpha_cop = m.asin((np.dot(Z23,r1))/(np.linalg.norm(Z23)*np.linalg.norm(r1)))

    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    r3_mag = np.linalg.norm(r3)

    # Make sure the vector have a moderate angular separation
    alpha12 = m.acos(np.dot(r1.flatten(),r2.flatten())/(r1_mag*r2_mag))
    alpha23 = m.acos(np.dot(r2.flatten(),r3.flatten())/(r2_mag*r3_mag))

    N = (r1_mag*Z23) + (r2_mag*Z31) + (r3_mag*Z12)
    D = Z12 + Z23 + Z31
    S = (r2_mag - r3_mag)*r1 + (r3_mag - r1_mag)*r2 + (r1_mag - r2_mag)*r3
    B = np.cross(D,r2.flatten())
    Lg = m.sqrt(mu/(np.linalg.norm(N)*np.linalg.norm(D)))
    v2 = (Lg/r2_mag)*B + (Lg*S.flatten())
    return v2,alpha_cop

def HerrickGibbs(r1,r2,r3,JD1,JD2,JD3):
    mu = 398600.4418
    delta_t31 = JD3 - JD1
    delta_t32 = JD3 - JD2
    delta_t21 = JD2 - JD1
    Z23 = np.cross(r2.flatten(),r3.flatten())

    # Must be low
    alpha_cop = (90*(np.pi/180)) - m.acos((np.dot(Z23,r1))/(np.linalg.norm(Z23)*np.linalg.norm(r1)))
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    r3_mag = np.linalg.norm(r3)

    # Make sure the vector have a moderate angular separation
    alpha12 = m.acos(np.dot(r1.flatten(),r2.flatten())/(r1_mag*r2_mag))
    alpha23 = m.acos(np.dot(r2.flatten(),r3.flatten())/(r2_mag*r3_mag))

    mu1 = (mu/(12*(r1_mag**3)))
    mu2 = (mu/(12*(r2_mag**3)))
    mu3 = (mu/(12*(r3_mag**3)))
    time1 = (1/(delta_t21*delta_t31))
    time2 = (1/(delta_t21*delta_t32))
    time3 = (1/(delta_t32*delta_t31))

    v2 = ((-delta_t32*(time1 + mu1)*r1) + 
            ((delta_t32 - delta_t21)*(time2 + mu2)*r2) +
            (delta_t21*(time3 + mu3)*r3))

    return v2,alpha_cop

def anglesOnlyGauss(L1,L2,L3,JD1,JD2,JD3,r_site1,r_site2,r_site3):
    mu = 398600.4418
    t1 = JD1 - JD2
    t3 = JD3 - JD2
    t1 = t1 * 86400
    t3 = t3 * 86400
    a1 = t3/(t3-t1)
    a1u = (t3*(((t3 - t1)**2) - (t3**2)))/(6*(t3 - t1))
    a3 = -t1/(t3-t1)
    a3u = (-t1*((((t3-t1)**2) - (t1**2))))/(6*(t3 - t1))
    L = np.array([L1, L2, L3])
    L = L.T.reshape((3,3))
    L_neg1 = np.array([L2[1]*L3[2] - L3[1]*L2[2], -L1[1]*L3[2] + L3[1]*L1[2], L1[1]*L2[2] - L2[1]*L1[2],
                        -L2[0]*L3[2] + L3[0]*L2[2], L1[0]*L3[2] - L3[0]*L1[2], -L1[0]*L2[2] + L2[0]*L1[2],
                        L2[0]*L3[1] - L3[0]*L2[1], -L1[0]*L3[1] + L3[0]*L1[1], L1[0]*L2[1] - L2[0]*L1[1]])/np.linalg.det(L)
    L_neg1 = L_neg1.reshape((3,3))
    L_neg1 = L_neg1.T

    r_site = np.array([r_site1, r_site2, r_site3])
    r_site = r_site.T.reshape((3,3))
    r_site1_mag = np.linalg.norm(r_site1)
    r_site2_mag = np.linalg.norm(r_site2)
    r_site3_mag = np.linalg.norm(r_site3)

    M = np.matmul(L_neg1,r_site)
    d1 = M[1][0]*a1 - M[1][1] + M[1][2]*a3
    d2 = M[1][0]*a1u + M[1][2]*a3u

    C = np.dot(L2.flatten(),r_site2.flatten())
    # Solve for roots for polynomial
    # r2 = symbols('r2')
    # equation = Eq((r2**8) - (((d1**2) + (2*C*d1) + (r_site2_mag**2))*(r2**6)) - (2*mu*((C*d2) + (d1*d2))*(r2**3)) - ((mu**2)*(d2**2)), 0)
    # solutions = solve(equation, r2)
    # valid_roots = [sol.evalf() for sol in solutions if (sol.is_real and (sol > 0))]
    # r2 = valid_roots[0]

    r2 = 10451.1381878519

    u = mu/(r2**3)
    c1 = a1 + a1u*u
    c2 = -1
    c3 = a3 + a3u*u
    Cnum = np.array([-c1, -c2, -c3])
    Cnum = Cnum.reshape((3,1))
    crho = np.matmul(M,Cnum)

    for ll in range(1,4):
        r1 = (crho[0]/c1)*L1 + r_site1
        r2 = (crho[1]/c2)*L2 + r_site2
        r3 = (crho[2]/c3)*L3 + r_site3
        r1_mag = np.linalg.norm(r1)
        r2_mag = np.linalg.norm(r2)
        r3_mag = np.linalg.norm(r3)
        # Below 1 deg Herrick-Gibbs is superior; above 5 deg Gibbs is superior
        v2,acop = Gibbs(r1,r2,r3)
        alpha12 = m.acos(np.dot(r1.flatten(),r2.flatten())/(r1_mag*r2_mag))
        alpha23 = m.acos(np.dot(r2.flatten(),r3.flatten())/(r2_mag*r3_mag))
        if acop < (np.pi/180):
            v2,acop = HerrickGibbs(r1,r2,r3,JD1 * 86400,JD2 * 86400,JD3 * 86400)
        (p,_,_,_,_,_,_) = RV2COE(r2,v2)

        if ll <= 2: # this section can be skipped
            u = mu/(r2_mag**3)
            rdot = np.dot(r2.flatten(),v2.flatten())/r2_mag
            udot = (-3*mu*rdot)/(r2_mag**4)
            tauSqrt = t1**2
            f1 = (1 - (0.5)*u*tauSqrt - (1/6)*udot*tauSqrt*t1
                    + (1.0/24.0)*u*u*tauSqrt*tauSqrt
                    + (1.0/30.0)*u*udot*tauSqrt*tauSqrt*t1)
            g1 = (t1 - (1.0/6.0)*u*t1*tauSqrt - (1.0/12.0)*udot*tauSqrt*tauSqrt
                    + (1.0/120.0)*u*u*tauSqrt*tauSqrt*t1
                    + (1.0/120.0)*u*udot*tauSqrt*tauSqrt*tauSqrt)
            tauSqrt = t3**2
            f3 = (1 - (0.5)*u*tauSqrt - (1/6)*udot*tauSqrt*t3
                    + (1.0/24.0)*u*u*tauSqrt*tauSqrt
                    + (1.0/30.0)*u*udot*tauSqrt*tauSqrt*t3)
            g3 = (t3 - (1.0/6.0)*u*t3*tauSqrt - (1.0/12.0)*udot*tauSqrt*tauSqrt
                    + (1.0/120.0)*u*u*tauSqrt*tauSqrt*t3
                    + (1.0/120.0)*u*udot*tauSqrt*tauSqrt*tauSqrt)
        else:
            f1 = 1 - (r1_mag/p)*(1 - m.cos(alpha12))
            f3 = 1 - (r3_mag/p)*(1 - m.cos(alpha23))
            g1 = (r1_mag*r2_mag*m.sin(-alpha12))/m.sqrt(mu*p)
            g3 = (r3_mag*r2_mag*m.sin(alpha23))/m.sqrt(mu*p)
        c1 = (g3)/((f1*g3) - (f3*g1))
        c3 = (-g1)/((f1*g3) - (f3*g1))
        Cnum_new = np.array([-c1, -c2, -c3])
        Cnum_new = Cnum_new.reshape((3,1))
        crho_new = np.matmul(M,Cnum_new)
        crho = crho_new
    r1 = (crho_new[0]/c1)*L1 + r_site1
    r2 = (crho_new[1]/c2)*L2 + r_site2
    r3 = (crho_new[2]/c3)*L3 + r_site3
    return r2,v2

def anglesDoubleGauss(L1,L2,L3,JD1,JD2,JD3,r_site1,r_site2,r_site3):
    mu = 398600.4418
    t1 = JD1 - JD2
    t3 = JD3 - JD2
    t1 = t1 * 86400
    t3 = t3 * 86400
    rsite1_mag = np.linalg.norm(r_site1)
    rsite2_mag = np.linalg.norm(r_site2)
    rsite3_mag = np.linalg.norm(r_site3)

    # Guess
    r1_mag = 12756.274
    r2_mag = 12820.05537
    c1 = 2*np.dot(L1.flatten(),r_site1.flatten())
    c2 = 2*np.dot(L2.flatten(),r_site2.flatten())
    for i in range(15):
        def F_Guess(r1_mag,r2_mag,c1,c2,t1,t3,mu):
            rho1 = (-c1 + m.sqrt((c1**2) - (4*((rsite1_mag**2) - (r1_mag**2)))))/2
            rho2 = (-c2 + m.sqrt((c2**2) - (4*((rsite2_mag**2) - (r2_mag**2)))))/2
            r1 = rho1*L1 + r_site1
            r2 = rho2*L2 + r_site2
            r1_mag = np.linalg.norm(r1)
            r2_mag = np.linalg.norm(r2)
            W = np.cross(r1.flatten(),r2.flatten())/(r1_mag*r2_mag)
            rho3 = (np.dot(-1*r_site3.flatten(),W))/(np.dot(L3.flatten(),W))
            r3 = rho3*L3 + r_site3
            r3_mag = np.linalg.norm(r3)

            delta_v21 = m.acos((np.dot(r2.flatten(),r1.flatten()))/(r2_mag*r1_mag))
            delta_v31 = m.acos((np.dot(r3.flatten(),r1.flatten()))/(r3_mag*r1_mag))
            delta_v32 = m.acos((np.dot(r3.flatten(),r2.flatten()))/(r3_mag*r2_mag))

            if (delta_v31 > np.pi):
                c1 = (r2_mag*m.sin(delta_v32))/(r1_mag*m.sin(delta_v31))
                c3 = (r2_mag*m.sin(delta_v21))/(r3_mag*m.sin(delta_v31))
                p = (c1*r1_mag + c3*r3_mag - r2_mag)/(c1 + c3 - 1)
            else:
                c1 = (r1_mag*m.sin(delta_v31))/(r2_mag*m.sin(delta_v32))
                c3 = (r1_mag*m.sin(delta_v21))/(r3_mag*m.sin(delta_v32))
                p = (c3*r3_mag - c1*r2_mag + r1_mag)/(-c1 + c3 + 1)
            ecosv1 = (p/r1_mag) - 1
            ecosv2 = (p/r2_mag) - 1
            ecosv3 = (p/r3_mag) - 1
            if (delta_v21 != np.pi):
                esinv2 = ((-m.cos(delta_v21)*ecosv2) + ecosv1)/m.sin(delta_v21)
            else:
                esinv2 = ((m.cos(delta_v32)*ecosv2) - ecosv3)/m.sin(delta_v31)
            e = m.sqrt((ecosv2**2) + (esinv2**2))
            a = p/(1-(e**2))
            n = m.sqrt(mu/a**3)
            S = (r2_mag/p)*m.sqrt(1 - (e**2))*esinv2
            C = (r2_mag/p)*((e**2) + ecosv2)
            delta_E32 = m.asin( ((r3_mag/m.sqrt(a*p))*m.sin(delta_v32)) - ((r3_mag/p)*(1-m.cos(delta_v32))*S) )
            delta_E21 = m.asin( ((r1_mag/m.sqrt(a*p))*m.sin(delta_v21)) + ((r1_mag/p)*(1-m.cos(delta_v21))*S) )
            delta_M32 = delta_E32 + (2*S*(m.sin(delta_E32/2)**2)) - (C*m.sin(delta_E32))
            delta_M12 = -delta_E21 + (2*S*(m.sin(delta_E21/2)**2)) + (C*m.sin(delta_E21))
            F1 = t1 - (delta_M12/n)
            F2 = t3 - (delta_M32/n)
            Q = m.sqrt((F1**2) + (F2**2))
            return F1,F2,a,delta_E32,r1,r2,r3

        # deltar1 = 0.005r1
        F1,F2,a,delta_E32,r1,r2,r3 = F_Guess(r1_mag,r2_mag,c1,c2,t1,t3,mu)
        F1_r1,F2_r1,_,_,_,_,_ = F_Guess(r1_mag + (0.005*r1_mag),r2_mag,c1,c2,t1,t3,mu)
        F1_r2,F2_r2,_,_,_,_,_ = F_Guess(r1_mag,r2_mag + (0.005*r2_mag),c1,c2,t1,t3,mu)
        deltaF1r1 = (F1_r1 - F1)/(0.005*r1_mag)
        deltaF2r1 = (F2_r1 - F2)/(0.005*r1_mag)
        deltaF1r2 = (F1_r2 - F1)/(0.005*r2_mag)
        deltaF2r2 = (F2_r2 - F2)/(0.005*r2_mag)

        delta = deltaF1r1*deltaF2r2 - deltaF2r1*deltaF1r2
        delta1 = deltaF2r2*F1 - deltaF1r2*F2
        delta2 = deltaF1r1*F2 - deltaF2r1*F1
        deltar1 = -delta1/delta
        deltar2 = -delta2/delta
        r1_mag = r1_mag + deltar1
        r2_mag = r2_mag + deltar2
    f = 1 - (a/r2_mag)*(1 - m.cos(delta_E32))
    g = t3 - m.sqrt((a**3)/mu)*(delta_E32 - m.sin(delta_E32))
    v2 = (r3 - f*r2)/g
    return r2,v2

# TODO: ADD Gooding's Method angles-only orbit determination

def lambertMinEnergy(r0,r):
    mu = 398600.4418
    r0_mag = np.linalg.norm(r0)
    r_mag = np.linalg.norm(r)
    cosdeltav = (np.dot(r0,r))/(r0_mag*r_mag)
    deltav = m.acos(cosdeltav)
    c = m.sqrt((r0_mag**2) + (r_mag**2) - (2*r0_mag*r_mag*cosdeltav))
    s = (r0_mag + r_mag + c)/2
    a_min = s/2
    p_min = (r0_mag*r_mag/c)*(1-cosdeltav)
    e_min = m.sqrt(1 - (2*p_min/s))
    alpha_e = np.pi
    beta_e = 2*m.asin(m.sqrt((s-c)/s))
    t_min_amin = m.sqrt((a_min**3)/mu)*(alpha_e - (beta_e - m.sin(beta_e)))
    t_min_abs = (1/3)*m.sqrt(2/mu)*((s**(3/2)) - ((s-c)**(3/2)))
    v0 = (m.sqrt(mu*p_min)/(r0_mag*r_mag*m.sin(deltav)))*(r - (1 - (r_mag/p_min)*(1 - cosdeltav))*r0)
    return a_min,e_min,t_min_abs,v0

def lambertGauss(r0,r,delta_t,tm):
    r0_mag = np.linalg.norm(r0)
    r_mag = np.linalg.norm(r)
    deltav = (np.dot(r0.flatten(),r.flatten()))
# _,_,_,_,_,_,_,_,_,JD_UT1_1,_ = convTime(2012,8,20,(11*3600) + (40*60) + 28.00,-0.609641,35.0)
# _,_,_,_,_,_,_,_,_,JD_UT1_2,_ = convTime(2012,8,20,(11*3600) + (48*60) + 28.00,-0.609641,35.0)
# _,_,_,_,_,_,_,_,_,JD_UT1_3,_ = convTime(2012,8,20,(11*3600) + (52*60) + 28.00,-0.609641,35.0)
# L1 = [m.cos(18.667717*(np.pi/180))*m.cos(0.939913*(np.pi/180)),m.cos(18.667717*(np.pi/180))*m.sin(0.939913*(np.pi/180)),m.sin(18.667717*(np.pi/180))]
# L1 = np.array(L1)
# L1 = L1.reshape((3,1))
# L2 = [m.cos(35.664741*(np.pi/180))*m.cos(45.025748*(np.pi/180)),m.cos(35.664741*(np.pi/180))*m.sin(45.025748*(np.pi/180)),m.sin(35.664741*(np.pi/180))]
# L2 = np.array(L2)
# L2 = L2.reshape((3,1))
# L3 = [m.cos(36.996583*(np.pi/180))*m.cos(67.886655*(np.pi/180)),m.cos(36.996583*(np.pi/180))*m.sin(67.886655*(np.pi/180)),m.sin(36.996583*(np.pi/180))]
# L3 = np.array(L3)
# L3 = L3.reshape((3,1))
# r_site1 = [4054.881,2748.195,4074.237]
# r_site1 = np.array(r_site1)
# r_site1 = r_site1.reshape((3,1))
# r_site2 = [3956.224,2888.232,4074.364]
# r_site2 = np.array(r_site2)
# r_site2 = r_site2.reshape((3,1))
# r_site3 = [3905.073,2956.935,4074.430]
# r_site3 = np.array(r_site3)
# r_site3 = r_site3.reshape((3,1))
# print(anglesDoubleGauss(L1,L2,L3,JD_UT1_1,JD_UT1_2,JD_UT1_3,r_site1,r_site2,r_site3))

deg2rad = np.pi/180
r = 6378.0

r_ITRF = np.array([-1033.4793830, 7901.2952754, 6380.3565958])
v_ITRF = np.array([-3.225636520, -2.872451450, 5.531924446])
r_ITRF = r_ITRF.reshape((3,1))
v_ITRF = v_ITRF.reshape((3,1))

#print(FK5(r_ITRF,v_ITRF,2004,4,6,(7*3600) + (51*60) + (28.386009),-0.4399619,32.0,-0.140682*(1/3600),0.333309*(1/3600),delta_delta_psi_1980=-0.052195*(1/3600),delta_delta_epsilon_1980=-0.003875*(1/3600),LOD=0.0015563))
#print(siteTrack(39.007*deg2rad,-104.883*deg2rad,2187/1000,604.68,205.6*deg2rad,30.7*deg2rad,2.08,0.15*deg2rad,0.17*deg2rad,1995,5,20,3*3600 + 17*60 +2, 0, 29, 0.0, 0.0))