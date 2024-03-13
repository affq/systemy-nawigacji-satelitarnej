import numpy as np
from alm_module import get_alm_data_str
import math

def julday(y,m,d,h=0):
    if m <= 2:
        y = y - 1
        m = m + 12
    jd = math.floor(365.25*(y+4716))+math.floor(30.6001*(m+1))+d+h/24-1537.5
    return jd

def get_gps_time(y,m,d,h=0,mnt=0,s=0):
    days = julday(y,m,d) - julday(1980,1,6)
    week = days//7
    day = days%7
    sow = day * 86400 + h * 3600 + mnt * 60 + s
    return week, sow

def flh2xyz(phi, lamb, h):
    '''
    funkcja przeliczająca współrzędne geograficzne na współrzędne kartezjańskie

    phi, lamb  w radianach
    '''
    N = a/np.sqrt(1-e2*np.sin(phi)**2)
    x = (N + h)*np.cos(phi)*np.cos(lamb)
    y = (N + h)*np.cos(phi)*np.sin(lamb)
    z = (N*(1-e2)+h)*np.sin(phi)
    return [x, y, z]

def Rneu(phi, lamb):
    '''
    funkcja zwracająca macierz rotacji z układu ECEF do układu NEU

    phi, lamb  w radianach
    '''
    R = np.array([[-np.sin(phi)*np.cos(lamb), -np.sin(lamb), np.cos(phi)*np.cos(lamb)],
                    [-np.sin(phi)*np.sin(lamb), np.cos(lamb), np.cos(phi)*np.sin(lamb)],
                    [np.cos(phi), 0, np.sin(phi)]])
    return R

a = 6378137 # wielka półoś elipsoidy GRS80 w metrach
e2 = 0.00669438002290 #kwadrat pierwszego mimośrodu dla elipsoidy GRS80

FI = 52 #szerokość geograficzna odbiornika w stopniach
LAMBDA = 21 #długość geograficzna odbiornika w stopniach
H = 100 #wysokość odbiornika w metrach

year, m, d = 2024, 2, 29 # czas na który chcemy wyznaczyć pozycję satelitów
h, mnt, s = 12, 0, 0

MI = 3.986005e14
OMEGA_E = 7.2921151467e-5

file = 'Almanac2024053.alm'

maska = 10 # maska elewacji w stopniach

satelites = []
A = []

WEEK, SOW = get_gps_time(year,m,d,h,mnt,s) # tydzień i czas w sekundach od początku tygodnia

nav, prn = get_alm_data_str(file)
satelity = nav[:,0]<100 # tylko giepeesy
nav = nav[satelity,:] 
prn = np.array(prn)[satelity]
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

odbiornik_x, odbiornik_y, odbiornik_z = flh2xyz(math.radians(FI), math.radians(LAMBDA), H) # współrzędne odbiornika w układzie xyz

print(f"Wyniki dla epoki y={year} m={m} d={d} h={h} mnt={mnt} s={s}")
print(f"Dane z pliku almanach: {file} \n")

print(f'Współrzędne odbiornika phi: {FI}, lambda: {LAMBDA} [deg]; h: {H}m')
print(f'Współrzędne XYZ odbiornika: {odbiornik_x}, {odbiornik_y}, {odbiornik_z}')
print(f'Maska elewacji: {maska} [deg]')

for sat in nav:
    x, y, z = satpos(SOW, WEEK, sat) # współrzędne satelity w układzie xyz
    sat_odb = np.array([x, y, z]) - np.array([odbiornik_x, odbiornik_y, odbiornik_z]) # wektor satelita-odbiornik xyz
    neu = Rneu(math.radians(FI), math.radians(LAMBDA)).T.dot(sat_odb) # wektor satelita-odbiornik neu
    az = math.degrees(math.atan2(neu[1], neu[0])) # azymut [deg]
    if az < 0:
        az += 360
    ro = np.sqrt(neu[0]**2 + neu[1]**2 + neu[2]**2) # odległość [m]
    el = math.degrees(math.asin(neu[2]/ro)) # elewacja [deg]

    print(f"SATELITA: {sat[0]}")
    print(f"Współrzędne satelity: {x}, {y}, {z} [m]")
    print(f"Wektor satelita-odbiornik w ukł. geocentrycznym Xsr: {sat_odb[0]}, {sat_odb[1]}, {sat_odb[2]} [m]")
    print(f"Wektor satelita-odbiornik w ukł. NEU: {neu[0]}, {neu[1]}, {neu[2]} [m]")
    print(f"Elewacja: {el} [deg]")
    print(f"Azymut: {az} [deg] \n")

    if el > maska:
        satelites.append(sat[0])
        A.append([-(x-odbiornik_x)/ro, -(y-odbiornik_y)/ro, -(z-odbiornik_z)/ro, 1])

A = np.array(A)
Q = np.linalg.inv(A.T.dot(A))
Qdiag = np.diag(Q)

GDOP = math.sqrt(Qdiag[0] + Qdiag[1] + Qdiag[2] + Qdiag[3])
PDOP = math.sqrt(Qdiag[0] + Qdiag[1] + Qdiag[2])
TDOP = math.sqrt(Qdiag[3])

Qneu = Rneu(math.radians(FI), math.radians(LAMBDA)).T.dot(Q[0:3,0:3]).dot(Rneu(math.radians(FI), math.radians(LAMBDA)))
Qneudiag = np.diag(Qneu)

HDOP = math.sqrt(Qneudiag[0] + Qneudiag[1])
VDOP = math.sqrt(Qneudiag[2])
PDOPneu = math.sqrt(Qneudiag[0] + Qneudiag[1] + Qneudiag[2])

print(f"GDOP: {round(GDOP, 5)}")
print(f"PDOP: {round(PDOP, 5)}")
print(f"TDOP: {round(TDOP, 5)}")
print(f"HDOP: {round(HDOP, 5)}")
print(f"VDOP: {round(VDOP, 5)}")

# o polnocy gop najnizszy, o 10 najwyzszy

