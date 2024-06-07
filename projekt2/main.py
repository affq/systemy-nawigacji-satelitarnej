import numpy as np
from readrnx_studenci import readrnxnav, readrnxobs, date2tow
from funcs import *
import math

OMEGA = 7.2921151467e-5 #[rad/s]
C = 299792458.0 #[m/s]

nav_file = 'rinex/BRDC00WRD_R_20240650000_01D_GN.rnx'
obs_file = 'rinex/JOZ200POL_R_20240650000_01D_30S_MO.rnx'

time_start = [2024, 3, 5, 0, 0, 0]  
time_end = [2024, 3, 5, 23, 59, 59]

obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
nav, inav = readrnxnav(nav_file)

zdrowe = nav[:,30] == 0
nav = nav[zdrowe,:]
inav = inav[zdrowe]

# współrzędne przybliżone odbiornika
xr0 = [3660000.,  1400000.,  5000000.]


el_mask = 10 # elevation mask/cut off in degrees
 

week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]

"""
Otwieramy dużą pętlę
for t in range(tow, tow_end+1, dt): gdzie dt równe 30
"""
t = 213300

index_t = iobs[:,2]==t
Pobs = obs[index_t,0]

# satelity = iobs[index_t,0]
satelity = np.array([25, 31, 32, 29, 28, 24, 20, 11, 12,  6])

tau = 0.07
dtr = 0
ro = 0

visible = False

for i in range(2):
    print(f"ITERACJA NUMER {i}")
    B, L, H = hirvonen(xr0[0], xr0[1], xr0[2])
    print(f"B, L, H dla iteracji {i} = {np.rad2deg(B), np.rad2deg(L), H}")
    A = []
    Y = []
    for nr, sat in enumerate(satelity):
        print(f"SATELITA G{sat}:")
        print(f"czas propagacji sygnału: {tau:.10f} s")
        idx_sat = inav == sat
        nav_sat = nav[idx_sat,:]
        dt = np.abs(t - nav_sat[:,17])
        idx_min_dt = np.argmin(dt)
        nav_choosen = nav_sat[idx_min_dt,:]

        ts = t - tau + dtr
        print(f"czas emisji sygnału: {ts} s")
        x0s, y0s, z0s, dt0s = satpos(ts, nav_choosen)
        print(f"współrzędne satelity na epokę ts: {x0s, y0s, z0s}")
        print(f"błąd zegara satelity na epokę ts: {dt0s} s")
        alfa = tau * OMEGA
        R = np.array([[np.cos(alfa), np.sin(alfa), 0], [-np.sin(alfa), np.cos(alfa), 0], [0, 0, 1]])
        print(f"macierz obrotu wokół osi Z: {R}")
        xyz0s = np.array([x0s, y0s, z0s])
        xs_rot = R@xyz0s
        print(f"współrzędne satelity po obrocie: {xs_rot}")
        sat_odb = np.array([xs_rot[0] - xr0[0], xs_rot[1] - xr0[1], xs_rot[2] - xr0[2]])
        neu = Rneu(B, L).T.dot(sat_odb)
        ro = np.sqrt((xs_rot[0] - xr0[0])**2 + (xs_rot[1] - xr0[1])**2 + (xs_rot[2] - xr0[2])**2)
        az = np.arctan2(neu[1], neu[0])
        el = np.arcsin(neu[2]/ro)
        print(f"odległość geometryczna ro: {ro}")
        print(f"wektor satelita-odbiornik: {sat_odb}")
        print(f"elewacja: {np.rad2deg(el)}°")
        
        # tau = ro/C
        
        if az < 0:
            az = az + 2*np.pi

        print(f"azymut: {np.rad2deg(az)}°")

        if el > np.radians(el_mask):
            visible = True
        else:
            visible = False

        """
                        Obliczamy poprawki atmosferyczne - dopiero wówczas, kiedy działać będzie nam program bez uwzględniania poprawek:
                            trop oraz iono
        """
        if visible:
            Pcalc = ro + c*dtr - c*dt0s
            print(f"Pcalc: {Pcalc}")
            y = Pobs[nr] - Pcalc
            print(f"element wektora wyrazów wolnych: {y} [m]")
            Y.append(y)
            a = [-(xs_rot[0] - xr0[0])/ro, -(xs_rot[1] - xr0[1])/ro, -(xs_rot[2] - xr0[2])/ro, 1]
            print(f"element macierzy A dla satelity: {a}")
            A.append(a)


    A = np.array(A)
    Y = np.array(Y)
    x = np.linalg.inv(A.T@A)@A.T@Y
    print(f"macierz A dla iteracji {i}: {A}")
    print(f"wektor wyrazów wolnych dla iteracji {i}: {Y}")
    print(f"wektor niewiadomych x dla iteracji {i}: {x}")

    xr = [xr0[0] + x[0], xr0[1] + x[1], xr0[2] + x[2]]
    print(f"współrzędne odbiornika dla iteracji {i}: \nx = {xr[0]}, \ny = {xr[1]}, \nz = {xr[2]}")

    dtr = dtr + x[3]/C
    print(f"błąd zegara odbiornika dla iteracji {i}: {dtr} s")
    exit()


"""
            Po skończeniu 5. iteracji, zbieramy obliczone współrzędne xr - warto zebrać również
            liczby obserwowanych satelitów, obliczone wartoci współczynników DOP (przynajmniej PDOP)
            
"""








