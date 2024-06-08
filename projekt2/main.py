import numpy as np
from readrnx_studenci import readrnxnav, readrnxobs, date2tow
from funcs import *

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
xr_reference = [3660000.,  1400000.,  5000000.]

el_mask = 10
 
week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]

satelity = []

tau = 0.07
dtr = 0
ro = 0
dt = 30

for t in range (tow, tow_end+1, dt):
    index_t = iobs[:,2]==t
    Pobs = obs[index_t,0]
    satelity = iobs[index_t,0]
    taus = []
    xrs_iter = []
    test_array = []
    for i in range(5):
        B, L, H = hirvonen(xr0[0], xr0[1], xr0[2])
        A = []
        Y = []
        for nr, sat in enumerate(satelity):
            if i == 0:
                tau = 0.07
            else:
                tau = taus[nr]
            idx_sat = inav == sat
            nav_sat = nav[idx_sat,:]
            dt = np.abs(t - nav_sat[:,17])
            idx_min_dt = np.argmin(dt)
            nav_choosen = nav_sat[idx_min_dt,:]
            ts = t - tau + dtr
            x0s, y0s, z0s, dt0s = satpos(ts, nav_choosen)
            alfa = tau * OMEGA
            R = np.array([[np.cos(alfa), np.sin(alfa), 0], [-np.sin(alfa), np.cos(alfa), 0], [0, 0, 1]])
            xyz0s = np.array([x0s, y0s, z0s])
            xs_rot = R@xyz0s
            sat_odb = np.array([xs_rot[0] - xr0[0], xs_rot[1] - xr0[1], xs_rot[2] - xr0[2]])
            neu = Rneu(B, L).T.dot(sat_odb)
            ro = np.sqrt((xs_rot[0] - xr0[0])**2 + (xs_rot[1] - xr0[1])**2 + (xs_rot[2] - xr0[2])**2)
            az = np.arctan2(neu[1], neu[0])
            el = np.arcsin(neu[2]/ro)
            taus.append(ro/C)
            if az < 0:
                az = az + 2*np.pi
            if el > np.radians(el_mask):
                tropo = Hopfield(H, np.rad2deg(el), i)
                iono = Klobuchar(np.rad2deg(B), np.rad2deg(L), np.rad2deg(el), np.rad2deg(az), ts)
                if i == 4:
                    test_array.append(f"sat = {sat} el = {np.rad2deg(el):.2f}° az = {np.rad2deg(az):.2f}° iono = {iono:.3f}m tropo = {tropo:.3f}m")
                # bez poprawek
                Pcalc = ro + c*dtr - c*dt0s 
                # z poprawkami
                # Pcalc = ro + c*dtr - c*dt0s + tropo + iono
                y = Pobs[nr] - Pcalc
                Y.append(y)
                a = [-(xs_rot[0] - xr0[0])/ro, -(xs_rot[1] - xr0[1])/ro, -(xs_rot[2] - xr0[2])/ro, 1]
                A.append(a)

        A = np.array(A)
        Y = np.array(Y)
        x = np.linalg.pinv(A.T@A)@A.T@Y
        xr0 = [xr0[0] + x[0], xr0[1] + x[1], xr0[2] + x[2]]
        dtr = dtr + x[3]/C
        xrs_iter.append(xr0)

    print(f"współrzędne odbiornika dla pierwszych 5 iteracji na czas {t}: {xrs_iter}")

for line in test_array:
    print(line)


"""
            Po skończeniu 5. iteracji, zbieramy obliczone współrzędne xr - warto zebrać również
            liczby obserwowanych satelitów, obliczone wartoci współczynników DOP (przynajmniej PDOP)
            
"""








