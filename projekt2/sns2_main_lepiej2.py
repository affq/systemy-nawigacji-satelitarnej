import numpy as np
from readrnx_studenci import readrnxnav, readrnxobs, date2tow
import math as mat

a = 6378137
e2 = 0.00669438002290

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

def Rneu(fi,lam):
    R = np.array(
            [
                [-np.sin(fi) * np.cos(lam), -np.sin(lam), np.cos(fi) * np.cos(lam)],
                [-np.sin(fi) * np.sin(lam), np.cos(lam), np.cos(fi) * np.sin(lam)],
                [np.cos(fi), 0, np.sin(fi)],

            ]
        )
    return R
def sat_pos(tow,nav_wybrane):
    u = 3.986005*10**14 #stała grawitacji 
    we = 7.2921151467*10**(-5) #prędkośc kątowa

    toe = nav_wybrane[17]
    alfa_f0 = nav_wybrane[6]
    alfa_f1 = nav_wybrane[7]
    alfa_f2 = nav_wybrane[8]
    gsp_week = nav_wybrane[27]
    a = nav_wybrane[16]**2
    e = nav_wybrane[14]
    i0 = nav_wybrane[21]
    omega0 = nav_wybrane[19]
    w = nav_wybrane[23]
    M0 = nav_wybrane[12]
    deltan = nav_wybrane[11]
    omega = nav_wybrane[24]
    idot = nav_wybrane[25]
    Cuc = nav_wybrane[13]
    Cus = nav_wybrane[15]
    Cic = nav_wybrane[18]
    Cis = nav_wybrane[20]
    Crc = nav_wybrane[22]
    Crs = nav_wybrane[10]

    # pkt1
    tk = tow-toe
    print('tk = ',tk)
    print('toe = ',toe)
    #print('tk = ',tk)

    n0 = np.sqrt(u / (a**3) )
    #('n0= ', n0)

    n = n0 + deltan
    #('n = ',n)

    Mk = M0 + n * tk
    #('Mk = ',Mk)

    epsilon = 1


    E = Mk
    while epsilon > 10**(-12):


        Enew = Mk + e* np.sin(E)
        epsilon = abs(E-Enew)
        E= Enew
        #('E = ',E)

    
    licznik = np.sqrt(1-e**2) * np.sin(E)
    mianownik = np.cos(E) - e
    vk = np.arctan2(licznik,mianownik)
    #('vk = ',vk)

    fi_k = vk + w
    #('fi _k= ',fi_k)
    
    delta_uk = Cus * np.sin(2*fi_k) + Cuc * np.cos(2*fi_k)
    delta_rk = Crs * np.sin(2*fi_k) + Crc * np.cos(2*fi_k)
    delta_ik = Cis * np.sin(2*fi_k) + Cic * np.cos(2*fi_k)
    #('delta_uk = ',delta_uk)
    #('delta_rk = ',delta_rk)
    #('delta_ik = ',delta_ik)

    uk = fi_k + delta_uk
    rk = a*(1-e*np.cos(E)) + delta_rk
    ik = i0 + idot * tk + delta_ik
    #('uk = ',uk)
    #('rk = ',rk)
    #('ik = ',ik)

    xk = rk * np.cos(uk)
    yk = rk * np.sin(uk)
    #('xk = ',xk)
    #('yk = ',yk)

    omega_k = omega0 + (omega - we) * tk - we * toe
    #('omega_k = ',omega_k)

    Xk = xk * np.cos(omega_k) - yk * np.cos(ik) * np.sin(omega_k)
    Yk = xk * np.sin(omega_k) + yk * np.cos(ik) * np.cos(omega_k)
    Zk = yk * np.sin(ik)
    #('Xk= ',Xk)
    #('Yk = ',Yk)
    #('Zk = ',Zk)

    dts = alfa_f0 + alfa_f1 * tk + alfa_f2 * tk**2
    #('dts = ',dts)

    c = 299792458
    c=c*c
    dtrel = (-2 * np.sqrt(u))/c * e * np.sqrt(a) * np.sin(E)
    dtsrel= dts + dtrel
    #('dtrel = ',dtrel)
    #('dtsrel = ',dtsrel)

    #nanarg

    XYZ = np.array([Xk, Yk, Zk])
    
    return XYZ, dtsrel

# cieżka do pliku nawigacyjnego
nav_file = 'rinex\BRDC00WRD_R_20240650000_01D_GN.rnx'
# cieżka do pliku obserwacyjnego
obs_file = 'rinex\JOZ200POL_R_20240650000_01D_30S_MO.rnx'

# zdefiniowanie czasu obserwacji: daty początkowej i końcowej
# dla pierwszej epoki z pliku będzie to:
time_start =  [2024, 3, 5, 0, 0, 0]  
time_end =    [2024, 3, 5, 23, 59, 59] 

# odczytanie danych z pliku obserwacyjnego
obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')     #mozna podac inny system np "E" - Galileo
# odczytanie danych z pliku nawigacyjnego:
nav, inav = readrnxnav(nav_file)

# iobs - nr satelity, sekunda od pocztku dnia, sekunda od pocztaku tygodnia
#obs - pseudoodleglosci



#usuwamy rekordy z chorymi satelitami
zdrowe = nav[:,30] == 0
nav = nav[zdrowe,:]
inav = inav[zdrowe]



#%%
"""
zdefiniowanie współrzędnych przybliżonych odbiornika - mogą to być współrzędne z nagłówka 
pliku obserwacyjnego, skopiowane "z palca" lub pobierane automatycznie z treci nagłówka pliku Rinex
"""
xr0 = [3660000.,  1400000.,  5000000.]

"""
Wprowadzenie ustawień, takich jak maska obserwacji, czy typ poprawki troposferycznej
"""
el_mask = 10 # elevation mask/cut off in degrees

"""
Przeliczenie daty początkowej i końcowej do sekund tygodnia GPS - niezbędne w celu
poprawnej definicji pętli związanej z czasem obserwacji w ciągu całej doby
"""
week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]
#%% Obliczenia

tow = 213300 #sekunda tygodnia dla godz 11.15 , potem zmienic z kazdym przejsciem petli
#1
index_t = iobs[:,2]==tow
Pobs = obs[index_t,0]
#2
satelity = np.array([25, 31, 32, 29, 28, 24, 20, 11, 12, 6])

#????????????????

#3
tau = 0.07
dtr = 0
tau_lista = []

for i in range(1):
    A=[]
    for sat in satelity:
        idx_sat = inav == sat
        nav_sat = nav[idx_sat,:]
        dt = np.abs(tow - nav_sat[:,17])
        idx_min_dt = np.argmin(dt)
        nav_choosen = nav_sat[idx_min_dt,:]
        
        
        ts = tow - tau + dtr
        #print ('ts = ',ts)
        xs0,dts = sat_pos(ts,nav_choosen)
        print ('współrzędne satelity xs0 = ',xs0)
        print ('dts = ',dts)

        #2. transformacja wsp satelity do ukladu chwilowego na moment odpioru sygnalu
        we = 7.2921151467*10**(-5)
        c = 299792458
        macierz = np.array([[np.cos(we*tau), np.sin(we*tau), 0],
                            [-np.sin(we*tau), np.cos(we*tau), 0],
                            [0, 0, 1]])
        xs = np.dot(macierz,xs0)
        #print ('współrzędne satelity w ukl. chwilowym xs = ',xs)

        #3. obliczenie pseudoodleglosi
        p0s = np.sqrt((xs[0]-xr0[0])**2 + (xs[1]-xr0[1])**2 + (xs[2]-xr0[2])**2)
        tau_new = p0s/c
        tau_lista.append(tau_new)
        #print('pseuoodleglosc = ',p0s)  #rho

        #wektor satelita odbiornik xsr
        xsr = xs - xr0
        #print('wektor satelita odbiornik xsr = ',xsr)

        #4. obliczenie kąta elewacji
        b,l,h = hirvonen(xr0[0],xr0[1],xr0[2])
        R = Rneu(b,l)
        xsr_neu = R.T.dot(xsr)
        el = np.arcsin(xsr_neu[2]/p0s)
        el_deg = el * 180 / np.pi
        #print('el = ',el_deg)
        az = np.arctan2(xsr_neu[1], xsr_neu[0])
        if az < 0:
            az = az + 2*np.pi
        az_deg = az * 180 / np.pi
        #print('az = ',az_deg)

        #pcalc
        pcalc = p0s - c*dts + c*dtr
        print('pcalc = ', pcalc)
        exit()