import numpy as np
from readrnx_studenci import readrnxnav, readrnxobs, date2tow
from funcs import *
import math

OMEGA = 7.2921151467e-5 #[rad/s]
C = 299792458.0 #[m/s]

# cieżka do pliku nawigacyjnego
nav_file = 'rinex/BRDC00WRD_R_20240650000_01D_GN.rnx'
# cieżka do pliku obserwacyjnego
obs_file = 'rinex/JOZ200POL_R_20240650000_01D_30S_MO.rnx'

# zdefiniowanie czasu obserwacji: daty początkowej i końcowej
# dla pierwszej epoki z pliku będzie to:
time_start = [2024, 3, 5, 0, 0, 0]  
time_end = [2024, 3, 5, 23, 59, 59]

# odczytanie danych z pliku obserwacyjnego
obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
# odczytanie danych z pliku nawigacyjnego:
nav, inav = readrnxnav(nav_file)


zdrowe = nav[:,30] == 0
nav = nav[zdrowe,:]
inav = inav[zdrowe]



# współrzędne przybliżone odbiornika
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

"""
Otwieramy dużą pętlę
for t in range(tow, tow_end+1, dt): gdzie dt równe 30
"""
t = 213300

print(week, t)

index_t = iobs[:,2]==t
Pobs = obs[index_t,0]
# satelity = iobs[index_t,0]
satelity = np.array([25, 31, 32, 29, 28, 24, 20, 11, 12,  6])

tau = 0.07
dtr = 0
ro = 0

A = []
Y = []
visible = False

for i in range(2):
    B, L, H = hirvonen(xr0[0], xr0[1], xr0[2])
    for sat in satelity:
        idx_sat = inav == sat
        nav_sat = nav[idx_sat,:]
        dt = np.abs(t - nav_sat[:,17])
        idx_min_dt = np.argmin(dt)
        nav_choosen = nav_sat[idx_min_dt,:]

        ts = t - tau + dtr
        x0s, y0s, z0s, dt0s = satpos(ts, nav_choosen)
        print(x0s, y0s, z0s, dt0s)
        alfa = tau * OMEGA
        R = np.array([[np.cos(alfa), np.sin(alfa), 0], [-np.sin(alfa), np.cos(alfa), 0], [0, 0, 1]])
        xyz0s = np.array([x0s, y0s, z0s])
        xs_rot = R@xyz0s
        sat_odb = np.array([xs_rot[0] - xr0[0], xs_rot[1] - xr0[1], xs_rot[2] - xr0[2]])
        neu = Rneu(B, L).T.dot(sat_odb)
        az = np.arctan2(neu[1], neu[0])
        ro = np.sqrt((xs_rot[0] - xr0[0])**2 + (xs_rot[1] - xr0[1])**2 + (xs_rot[2] - xr0[2])**2)
        tau = ro/C
        
        if az < 0:
            az = az + 2*np.pi
        el = np.arcsin(neu[2]/ro)

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
            y = Pobs - Pcalc
            Y.append(y)
            a = [-(xs_rot[0] - xr0[0])/ro, -(xs_rot[1] - xr0[1])/ro, -(xs_rot[2] - xr0[2])/ro, 1]
            A.append(a)
        exit()

A = np.array(A)
Y = np.array(Y)
x = np.linalg.lstsq(A, Y, rcond=None)[0]

xr = [xr0[0] + x[0], xr0[1] + x[1], xr0[2] + x[2]]
dtr = dtr + x[3]/C

"""
            Po skończeniu 5. iteracji, zbieramy obliczone współrzędne xr - warto zebrać również
            liczby obserwowanych satelitów, obliczone wartoci współczynników DOP (przynajmniej PDOP)
            
"""








