import numpy as np
from alm_module import get_alm_data_str
import math

def julday(y,m,d,h=0):
    if m <= 2:
        y = y - 1
        m = m + 12
    jd = math.floor(365.25*(y+4716))+math.floor(30.6001*(m+1))+d+h/24-1537.5;
    return jd

def get_gps_time(y,m,d,h=0,mnt=0,s=0):
    days = julday(y,m,d) - julday(1980,1,6)
    week = days//7
    day = days%7
    sow = day * 86400 + h * 3600 + mnt * 60 + s
    return week, sow

a = 6378137 # wielka półoś elipsoidy GRS80 w metrach
e2 = 0.00669438002290 # kwadrat pierwszego mimośrodu dla elipsoidy GRS80

def flh2xyz(phi, lamb, h):
    N = a/np.sqrt(1-e2*np.sin(phi)**2)
    x = (N + h)*np.cos(phi)*np.cos(lamb)
    y = (N + h)*np.cos(phi)*np.sin(lamb)
    z = (N*(1-e2)+h)*np.sin(phi)
    return [x, y, z]

def Rneu(phi, lamb):
    R = np.array([[-np.sin(phi)*np.cos(lamb), -np.sin(lamb), np.cos(phi)*np.cos(lamb)],
                    [-np.sin(phi)*np.sin(lamb), np.cos(lamb), np.cos(phi)*np.sin(lamb)],
                    [np.cos(phi), 0, np.sin(phi)]])
    return R

FI = 52
LAMBDA = 21
H = 100

FI, LAMBDA, H = flh2xyz(FI, LAMBDA, H)

nav, prn = get_alm_data_str('Almanac2024053.alm')

satelity = nav[:,0]<400
nav = nav[satelity,:]
prn = np.array(prn)[satelity]

y, m, d = 2024, 2, 29
h, mnt, s = 12, 0, 0

week, sow = get_gps_time(y,m,d,h,mnt,s)
print('Tydzień:', week, 'Sekunda tygodnia:', sow)

MI = 3.986005e14 
OMEGA_E = 7.2921151467e-5

wiersz_nav = nav[0,:]

def satpos(sow, week, wiersz_nav):
    e = wiersz_nav[2] # ekscentryczność orbity (mimośród)
    a = wiersz_nav[3]**2 # półoś wielka [m]
    Omega = wiersz_nav[4] * np.pi/180 # długość węzła wstępującego [rad]
    omega = wiersz_nav[5] * np.pi/180 # argument perigeum [rad]
    M0 = wiersz_nav[6] * np.pi/180 # anomalia średnia [rad]
    ToA = wiersz_nav[7] # czas odniesienia almanachu [s]
    i = (54 + wiersz_nav[8]) * np.pi/180 # inklinacja [rad]
    Omega_dot = (wiersz_nav[9]/1000) * np.pi/180 # tempo zmiany rektanscenzji węzła wstępującego [rad/s]
    gps_week = wiersz_nav[12]

    t_z_tygodniami = week * 7 * 86400 + sow
    ToA_z_tygodniami = gps_week * 7 * 86400 + ToA

    tk = t_z_tygodniami - ToA_z_tygodniami

    n = math.sqrt(MI/a**3) # Wyznaczenie średniej prędkości kątowej n, znanej jako ruch średni (ang. mean motion) na podstawie III prawa Kepplera

    Mk = M0 + n * tk # Poprawiona anomalia średnia na epokę tk

    #Wyznaczenie anomalii mimośrodowej (Równanie Kepplera)
    error = 1
    E = Mk
    while error > 1e-12:
        Enew = Mk + e * math.sin(E)
        error = abs(Enew - E)
        E = Enew
    
    #Wyznaczenie anomalii prawdziwej
    vk = math.atan2(math.sqrt(1-e**2)*math.sin(E), math.cos(E)-e)
    
    #Wyznaczenie argumentu szerokości
    Phik = vk + omega

    #Wyznaczenie promienia orbity
    rk = a * (1 - e * math.cos(E))

    #Wyznaczenie pozycji satelity w układzie orbity
    xk = rk * math.cos(Phik)
    yk = rk * math.sin(Phik)

    #Poprawiona długość węzła wstępującego
    Omega_k = Omega + (Omega_dot - OMEGA_E) * tk - (OMEGA_E * ToA)

    #Wyznaczenie pozycji satelity w układzie geocentrycznym ECEF:
    x = xk * math.cos(Omega_k) - yk * math.cos(i) * math.sin(Omega_k)
    y = xk * math.sin(Omega_k) + yk * math.cos(i) * math.cos(Omega_k)
    z = yk * math.sin(i)

    return x, y, z

positions = []
maska = 10

for i in range(0, 24):
    for sat in nav:
        x, y, z = satpos(sow, week, sat)
        xryrzr = np.array([x, y, z]) - np.array([FI, LAMBDA, H])
        neu = Rneu(FI, LAMBDA).dot(xryrzr)
        az = math.atan2(neu[0], neu[1])
        el = math.asin(neu[2]/np.sqrt(neu[0]**2 + neu[1]**2 + neu[2]**2))

        if el > maska:
            #wiersz A...
            A = []

A = np.array(A)
Q = np.linalg.inv(A.T.dot(A))
Q = np.diag(Q)
GDOP = math.sqrt(np.trace(Q))
PDOP = math.sqrt(Q[0] + Q[1] + Q[2])
HDOP = math.sqrt(Q[0] + Q[1])
VDOP = math.sqrt(Q[2])

Qneu = Rneu(FI, LAMBDA).T.dot(Q[0:3]).dot(Rneu(FI, LAMBDA))

print(positions[0])
print(positions[1])
print(positions[2])




# o polnocy gop najnizszy, o 10 najwyzszy

# dla satelitow > maska liczymy wiersz macierzy A