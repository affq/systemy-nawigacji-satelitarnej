import wget
import numpy as np
from alm_module import get_alm_data_str
import math

# file = 'ftp://ftp.trimble.com/pub/eph/Almanac.alm'
# wget.download(file, 'Almanac.alm')

nav, prn = get_alm_data_str('Almanac2024053.alm')

satelity = nav[:,0]<400
nav = nav[satelity,:]
prn = np.array(prn)[satelity]

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

y, m, d = 2024, 2, 29
h, mnt, s = 12, 0, 0

week, sow = get_gps_time(y,m,d)
print('Tydzień:', week, 'Sekunda tygodnia:', sow)

MI = 3.986005e14 
OMEGA_E = 7.2921151467e-5

wiersz_nav = nav[0,:]

def satpos(t, week, wiersz_nav):
    e = wiersz_nav[2] # ekscentryczność
    a = wiersz_nav[3]**2 # półoś wielka [m]
    Omega = wiersz_nav[4] * np.pi/180 # długość węzła wstępującego [rad]
    omega = wiersz_nav[5] * np.pi/180 # argument perigeum [rad]
    M0 = wiersz_nav[6] * np.pi/180 # anomalia średnia [rad]
    ToA = wiersz_nav[7] # czas almanac [s]
    i = (54 + wiersz_nav[8]) * np.pi/180 # inklinacja [rad]


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
    

    return x, y, z