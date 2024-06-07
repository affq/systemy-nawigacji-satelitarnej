import math as mat
import numpy as np

a = 6378137
e2 = 0.00669438002290
mi = 3.986005e14
omegae = 7.2921151467e-5
c = 299792458.0

def Np(B):
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N

def hirvonen(X,Y,Z):
    r = (X**2 + Y**2)**0.5
    B = mat.atan(Z/(r*(1-e2)))
    
    while 1:
        N = Np(B)
        H = r/np.cos(B) - N
        Bst = B
        B = mat.atan(Z/(r*(1-(e2*(N/(N+H))))))    
        if abs(Bst-B)<(0.00001/206265):
            break
    L = mat.atan2(Y,X)
    N = Np(B)
    H = r/np.cos(B) - N
    return B, L, H 

def Rneu(phi, lamb):
    R = np.array([[-np.sin(phi)*np.cos(lamb), -np.sin(lamb), np.cos(phi)*np.cos(lamb)],
                    [-np.sin(phi)*np.sin(lamb), np.cos(lamb), np.cos(phi)*np.sin(lamb)],
                    [np.cos(phi), 0, np.sin(phi)]])
    return R

def satpos(t, nav):
    gps_week = nav[27]
    toe = nav[17]
    alfaf0 = nav[6] # współczynnik wielomianu do poprawki zegara satelity (opóźnienie) [s]
    alfaf1 = nav[7] # współczynnik wielomianu do poprawki zegara (dryft) [s/s]
    alfaf2 = nav[8] # współczynnik wielomianu do poprawki zegara (częstotliwość dryftowania) [s/s^2]

    #elementy orbity keplerowskiej
    a = nav[16]**2
    e = nav[14]
    i0 = nav[21]
    Omega0 = nav[19]
    omega = nav[23]
    M0 = nav[12]

    # parametry perturbacyjne
    Deltan = nav[11]
    Omegadot = nav[24]
    i0dot = nav[25]
    Cuc = nav[13]
    Cus = nav[15]
    Cic = nav[18]
    Cis = nav[20]
    Crc = nav[22]
    Crs = nav[10]

    # czas, jaki upłynął od epoki wyznaczenia efemerydy
    tk = t - toe
    print(f"t = {t}")
    print(f"toe = {toe}")
    print(f"tk = {tk}")

    # średnia prędkość kątowa (ruch średni)
    n0 = (mi/a**3)**0.5

    # poprawiony ruch średni
    n = n0 + Deltan

    # poprawiona anomalia średnia na epokę tk
    Mk = M0 + n * tk

    # anomalia mimośrodowa
    error = 1
    Ek = Mk
    while error > 1e-12:
        E_new = Mk + e * np.sin(Ek)
        error = abs(E_new - Ek)
        Ek = E_new

    # anomalia prawdziwa
    vk = np.arctan2((1 - e**2)**0.5 * np.sin(Ek), np.cos(Ek) - e)

    # argument szerokości
    Phik = vk + omega

    # poprawka do argumentu szerokości
    Deltauk = Cus * np.sin(2 * Phik) + Cuc * np.cos(2 * Phik)

    # poprawka do promienia orbity
    Deltark = Crs * np.sin(2 * Phik) + Crc * np.cos(2 * Phik)

    # poprawka do inklinacji orbity
    Deltaik = Cis * np.sin(2 * Phik) + Cic * np.cos(2 * Phik)

    # poprawiony argument szerokości uk
    uk = Phik + Deltauk

    # poprawiony promień orbity rk
    rk = a * (1 - e * np.cos(Ek)) + Deltark

    # poprawiona inklinacja orbity ik
    ik = i0 + i0dot * tk + Deltaik

    # pozycja satelity w układzie orbity
    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)

    if abs(rk - (xk**2 + yk**2)**0.5) > 1e-2:
        raise ValueError("coś nie tak w układzie orbity")

    # poprawiona długość węzła wstępującego
    Omegak = Omega0 + (Omegadot - omegae) * tk - omegae * toe

    # pozycja satelity w układzie ecef
    x = xk * np.cos(Omegak) - (yk * np.cos(ik) * np.sin(Omegak))
    y = xk * np.sin(Omegak) + (yk * np.cos(ik) * np.cos(Omegak))
    z = yk * np.sin(ik)

    if abs((x**2 + y**2 + z**2)**0.5 - rk) > 1e-2:
        raise ValueError("coś nie tak w układzie ecef")

    # błąd synchronizacji zegara satelity
    deltats = alfaf0 + alfaf1 * (t - toe) + alfaf2 * (t - toe)**2

    # poprawka relatywistyczna
    deltatrel = (-2 * (mi)**0.5 / c**2) * e * a**0.5 * np.sin(Ek)

    # błąd synchronizacji zegara satelity z poprawką relatywistyczną
    deltatsrel = deltats + deltatrel

    return x, y, z, deltatsrel

def Hopfield(h, el, mask, iteration):
    if el > mask:
        if iteration==0:
            dT = 0
    else:
        hort = h - 31.36
        p = 1013.25 * (1 - 0.0000226*hort)**5.225
        temp = 291.15 - 0.0065 * hort
        Rh = 0.5 * np.exp(-0.0006396*hort)
        e = 6.11 * Rh * 10**((7.5*(temp-273.15))/(temp - 35.85))

        Nd0 = 77.64*p/temp
        Nw0 = -12.96*e/temp + 3.718*10**5*e/temp**2
        hd = 40136 + 148.72 * (temp - 273.15)
        hw = 11000
        dTd0 = 10**(-6)/5 * Nd0 * hd
        dTw0 = 10**(-6)/5 * Nw0 * hw

        md = 1/(np.sin(np.deg2rad(np.sqrt(el**2 + 6.25))))
        mw = 1/(np.sin(np.deg2rad(np.sqrt(el**2 + 2.25))))
        dT = dTd0*md + dTw0*mw 

alfa = [2.4214E-08, 0.0000E+00, -1.1921E-07, 5.9605E-08]
beta = [1.2902E+05, 0.0000E+00, -1.9661E+05, -6.5536E+04]

def Klobuchar(phi, lamb, el, az, tgps, alfa, beta):
    phi = 52
    lamb = 21
    el = 30
    az = 180
    tgps = 43200

    phis = phi/180
    lambs = lamb/180
    els = el/180
    azs = az/180


    # 1. kąt geocentryczny
    psi = 0.0137 / (els + 0.11) - 0.022

    # 2. szerokość geograficzna IPP
    phi_ipp = phis  + psi * np.cos(np.deg2rad(az))

    if phi_ipp > 0.416:
        phi_ipp = 0.416
    elif phi_ipp < -0.416:
        phi_ipp = -0.416

    # 3. długość geograficzna IPP
    lamb_ipp = lambs + (psi * np.sin(np.deg2rad(az)) / np.cos(phi_ipp*np.pi))

    # 4. szerokość geomagnetyczna IPP
    phim = phi_ipp + 0.064 * np.cos((lamb_ipp - 1.617) * np.pi)