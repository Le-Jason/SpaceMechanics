import math as m
import numpy as np

def c2c3(psi):
    if psi > 0.000001:
        c2 = (1 - m.cos(m.sqrt(psi)))/psi
        c3 = (m.sqrt(psi) - m.sin(m.sqrt(psi)))/(m.sqrt(psi**3))
    else:
        if (psi < -0.000001):
            c2 = (1 - m.cosh(m.sqrt(-psi)))/psi
            c3 = (m.sinh(m.sqrt(-psi)) - m.sqrt(-psi))/(m.sqrt((-psi)**3))
        else:
            c2 = 1/2
            c3 = 1/6
    return (c2,c3)

def KepEqtnE(M,e):
    tol = 0.00000001
    if (M > m.pi) or ((-m.pi < M) and (M < 0)):
        En = M - e
    else:
        En = M + e
    E = En + (M - En + (e*m.sin(En)))/(1 - (e*m.cos(En)))
    while (abs(E - En) > tol):
        En = E
        E = En + (M - En + (e*m.sin(En)))/(1 - (e*m.cos(En)))
    return E

def cot(x):
    return 1 / m.tan(x)

def arccot(x):
    return (m.pi / 2) - m.atan(x)

def KepEqtnP(deltat,p):
    mu = 398600.4418
    np = 2 * m.sqrt(mu/(p**3))
    cot2s = (3/2)*np*deltat
    s = arccot(cot2s)/2
    tan3w = m.tan(s)
    w = m.atan(tan3w**(1./3.))
    B = 2*cot(2*w)
    return B

def KepEqtnH(M,e):
    tol = 0.00000001
    if (e < 1.6):
        if ((M > m.pi) or ((-m.pi < M) and (M < 0))):
            Hn = M - e
        else:
            Hn = M + e
    else:
        if ((e < 3.6) and (abs(M) > m.pi)):
            Hn = M - (m.copysign(1,M)*e)
        else:
            Hn = M / (e - 1)
    H = Hn + (M - (e*m.sinh(Hn)) + Hn)/((e*m.cosh(Hn)) - 1)
    while (abs(H - Hn) > tol):
        Hn = H
        H = Hn + (M - (e*m.sinh(Hn)) + Hn)/((e*m.cosh(Hn)) - 1)
    return H

def v2Anomaly(e,v):
    if (e < 1.0):
        # sinE = (m.sin(v)*m.sqrt(1 - (e**2)))/(1 + (e*m.cos(v)))
        cosE = (e + m.cos(v)) / (1 + (e*m.cos(v)))
        E = m.acos(cosE)
        return E
    elif ( e == 1.0):
        B = m.tan(v/2)
        return B
    elif (e > 1.0):
        # sinhH = (m.sin(v)*m.sqrt((e**2) - 1))/(1 + (e*m.cos(v)))
        coshH = (e + m.cos(v)) / (1 + (e*m.cos(v)))
        H = m.acosh(coshH)
        return H
    print("v2Anomaly : error")
    return 1

def Anomaly2v(e,E,p=0,r=0):
    if (e < 1.0):
        # sinv = (m.sin(E)*m.sqrt(1 - (e**2)))/(1 - (e*m.cos(v)))
        cosv = (m.cos(E) - e) / (1 - (e*m.cos(E)))
        v = m.acos(cosv)
        return v
    elif ( e == 1.0):
        # sinv = (p*E)/r
        cosv = (p - r)/r
        v = m.acos(cosv)
        return v
    elif (e > 1.0):
        # sinhH = (-m.sinh(E)*m.sqrt((e**2) - 1))/(1 - (e*m.cosh(H)))
        coshv = (m.cosh(E) - e) / (1 - (e*m.cosh(E)))
        v = m.cos(coshv)
        return v
    print("Anomaly2v : error")
    return 1

def RV2COE(r,v):
    mu = 398600.4418
    h = np.cross(r,v)
    hNorm = m.sqrt(np.dot(h,h))
    K = np.array([0,0,1])
    n = np.cross(K,h)
    nNorm = m.sqrt(np.dot(n,n))
    vNorm = m.sqrt(np.dot(v,v))
    rNorm = m.sqrt(np.dot(r,r))
    e = ((((vNorm**2) - (mu/rNorm))*r) - (np.dot(r,v)*v)) / mu
    l = ((vNorm**2)/2) - (mu/rNorm)
    eNorm = m.sqrt(np.dot(e,e))
    if (eNorm != 1.0):
        a = -mu/(2*l)
        p = a*(1 - (eNorm**2))
    else:
        p = (hNorm**2)/mu
        a = 10000000000000000000 #inf
    i = m.acos(h[2]/hNorm)
    if (n[1] >= 0):
        RAAN = m.acos(n[0]/nNorm)
    elif (n[1] < 0):
        RAAN = (2*m.pi) - m.acos(n[0]/nNorm)
    if (e[2] >= 0):
        w = m.acos(np.dot(n,e)/(nNorm*eNorm))
    elif (e[2] < 0):
        w = (2*np.pi) - m.acos(np.dot(n,e)/(nNorm*eNorm))
    if (np.dot(r,v) >= 0):
        truev = m.acos(np.dot(e,r)/(eNorm*rNorm))
    elif (np.dot(r,v) < 0):
        truev = (2*np.pi) - m.acos(np.dot(e,r)/(eNorm*rNorm))
    if (e[1] >= 0):
        wtrue = m.acos(e[0]/eNorm)
    elif (e[1] < 0):
        wtrue = (2*np.pi) - m.acos(e[0]/eNorm)
    if (r[2] >= 0):
        u = m.acos(np.dot(n,r)/(nNorm*rNorm))
    elif (r[2] < 0):
        u = (2*np.pi) - m.acos(np.dot(n,r)/(nNorm*rNorm))
    if (r[1] >= 0):
        lambdatrue = m.acos(r[0]/rNorm)
    elif (r[1] < 0):
        lambdatrue = (2*np.pi) - m.acos(r[0]/rNorm)
    return (p,a,eNorm,i,RAAN,w,truev)

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



def COE2RV(a,e,i,RAAN,w,v,u=0,lambdatrue=0,wtrue=0):
    mu = 398600.4418
    p = a*(1 - e**2)
    # if circular equatorial
    # w = 0.0
    # RAAN = 0.0
    # v = lambdatrue

    # if circular inclined
    # w = 0.0
    # v = u

    # if elliptical equatorial
    # RAAN = 0.0
    # w = wtrue

    rPQW = [(p*m.cos(v))/(1 + (e*m.cos(v))),
            (p*m.sin(v))/(1 + (e*m.cos(v))),
            0.0]
    vPQW = [-m.sqrt(mu/p)*m.sin(v),
            m.sqrt(mu/p)*(e + m.cos(v)),
            0.0]
    
    # 3-1-3 Rotation ROT3(-RAAN)*ROT1(-i)*ROT3(-w)
    T_PQW_IJK = np.matmul(ROT1(-i), ROT3(-w))
    T_PQW_IJK = np.matmul(ROT3(-RAAN), T_PQW_IJK)

    rIJK = np.matmul(T_PQW_IJK,rPQW)
    vIJK = np.matmul(T_PQW_IJK,vPQW)

    return (rIJK,vIJK)

def KeplerCOE(r0,v0,deltat):
    [p,a,e,i,RAAN,w,truev] = RV2COE(r0,v0)
    mu = 398600.4418
    np = m.sqrt(mu/a**3)
    if (e != 0):
        E0 = v2Anomaly(e,truev)
    else:
        E0 = truev # E = u or E = lambdatrue
    if (e < 1.0):
        M0 = E0 - (e*m.sin(E0))
        M = M0 + (np*deltat)
        E = KepEqtnE(M,e)
    elif (e == 1.0):
        h = np.corss(r0,v0)
        hNorm = m.sqrt(np.dot(h,h))
        p = (hNorm**2) / mu
        M0 = E0 + (E0**3/3)
        E = KepEqtnP(deltat,p)
    elif (e > 1.0):
        M0 = (e*m.sinh(E0) - E0)
        M = M0 + (np*deltat)
        E = KepEqtnH(M,e)
    if (e == 0):
        cur_v = Anomaly2v(e,E,p,m.sqrt(np.dot(r0,r0)))
    elif (e != 0):
        cur_v = Anomaly2v(e,E)
    else:
        u = E
        lambdatrue = E
    [r,v] = COE2RV(a,e,i,RAAN,w,cur_v)
    return (r,v)

def Kepler(r0,v0,deltat):
    mu = 398600.4418
    v0Norm = m.sqrt(np.dot(v0,v0))
    r0Norm = m.sqrt(np.dot(r0,r0))
    l = ((v0Norm**2)/2) - (mu/r0Norm)
    alpha = (-(v0Norm**2)/mu) + (2/r0Norm)
    if (alpha > 0.000001):
        Xn = m.sqrt(mu) * deltat * alpha
    elif (alpha == 1.0):
        print("Kepler : Too close to converge")
        return 0
    elif (abs(alpha) < 0.000001):
        h = np.cross(r0,v0)
        hNorm = m.sqrt(np.dot(h,h))
        p = (hNorm**2) / mu
        n = 2 * m.sqrt(mu/(p**3))
        cot2s = (3/2)*n*deltat
        s = arccot(cot2s)/2
        tan3w = m.tan(s)
        w = m.atan(tan3w**(1./3.))
        Xn = m.sqrt(p)*2*cot(2*w)
    elif (alpha < -0.000001):
        a = 1/alpha
        top = -2*mu*alpha*deltat
        bot = np.dot(r0,v0) + (m.copysign(1,deltat)*m.sqrt(-mu*a)*(1 - (r0Norm*alpha)))
        Xn = m.copysign(1,deltat)*m.sqrt(-a)*np.log(top/bot)
    psi = (Xn**2) * alpha
    [c2,c3] = c2c3(psi)
    rtemp = ((Xn**2)*c2) + (np.dot(r0,v0)/m.sqrt(mu))*Xn*(1-(psi*c3)) + (r0Norm*(1 - (psi*c2)))
    X = Xn + ((m.sqrt(mu)*deltat) - ((Xn**3)*c3) - ((np.dot(r0,v0)/m.sqrt(mu))*(Xn**2)*c2) - (r0Norm * Xn * (1 - (psi*c3)))) / rtemp
    tol = 0.000001
    while (abs(X - Xn) > tol):
        Xn = X
        psi = (Xn**2) * alpha
        [c2,c3] = c2c3(psi)
        rtemp = ((Xn**2)*c2) + (np.dot(r0,v0)/m.sqrt(mu))*Xn*(1-(psi*c3)) + (r0Norm*(1 - (psi*c2)))
        X = Xn + ((m.sqrt(mu)*deltat) - ((Xn**3)*c3) - ((np.dot(r0,v0)/m.sqrt(mu))*(Xn**2)*c2) - (r0Norm * Xn * (1 - (psi*c3)))) / rtemp
    f = 1 - ((X**2)/r0Norm)*c2
    g = deltat - ((X**3)/m.sqrt(mu))*c3
    fdot = (m.sqrt(mu)/(rtemp*r0Norm))*X*((psi * c3) - 1)
    gdot = 1 - (((X**2)/rtemp)*c2)
    r = (f*r0) + (g*v0)
    v = (fdot*r0) + (gdot*v0)
    if (m.isclose(((f*gdot) - (fdot*g)), 1, abs_tol = 1e-8)):
        return r,v
    else:
        print("Kepler : fgdot - fdotg did not equal 1")
        return r,v
    
def findTOF(r0,r,p):
    mu = 398600.4418
    r0Norm = m.sqrt(np.dot(r0,r0))
    rNorm = m.sqrt(np.dot(r,r))
    cosdeltav = (np.dot(r0,r))/(r0Norm * rNorm)
    deltav = m.acos(cosdeltav)
    k = r0Norm*rNorm*(1 - cosdeltav)
    l = r0Norm + rNorm
    ml = r0Norm*rNorm*(1 + cosdeltav)
    a = (ml*k*p)/((((2*ml) - (l**2))*p**2) + (2*k*l*p) - (k**2))
    f = 1 - ((rNorm/p)*(1 - cosdeltav))
    g = (r0Norm * rNorm * m.sin(deltav))/(m.sqrt(mu*p))
    if (a > 0.0):
        fdot = m.sqrt(mu/p)*m.tan(deltav/2)*(((1-cosdeltav)/p) - (1/r0Norm) - (1/rNorm))
        cosdeltaE = 1 - ((r0Norm/a)*(1 - f))
        sindeltaE = (-r0Norm*rNorm*fdot)/(m.sqrt(mu*a))
        deltaE = m.acos(cosdeltaE)
        TOF = g + (m.sqrt(a**3/mu)*(deltaE - sindeltaE))
    if (a > 100000000000000000000):
        c = m.sqrt((r0Norm**2) + (rNorm**2) - (2*r0Norm*rNorm*m.cos(deltav)))
        s = (r0Norm + rNorm + c) / 2
        TOF = (2/3)*(np.sqrt(s**3/(2*mu)))*(1 - (((s-c)/s)**(3/2)))
    elif (a < 0.0):
        coshdeltaH = 1 + ((f - 1)*(r0Norm/a))
        deltaH = m.acosh(coshdeltaH)
        TOF = g + ((m.sqrt(((abs(a))**3)/mu))*(m.sinh(deltaH) - deltaH))
    return TOF

# deg2rad = m.pi/180
# print(findTOF(np.array([-2574.9533,4267.0671,4431.5025]),np.array([2700.6738,-4303.5378,-4358.2499]),6681.571))