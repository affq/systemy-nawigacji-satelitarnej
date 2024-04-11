from readrnx_studenci import readrnxnav, date2tow
import numpy as np

nav_file = r"BRDC00WRD_R_20240650000_01D_GN.rnx"
nav, inav = readrnxnav(nav_file)
zdrowe = nav[:, 30] == 0
nav = nav[zdrowe,:]
inav = inav[zdrowe]

def satpos_fromrinex():
    pass

# wejście
data_obliczen = [2024, 3, 5, 11, 15, 0]
week, tow, dow = date2tow(data_obliczen)

sat = 2

indeks_satelity = inav == sat
nav_satelity = nav[indeks_satelity, :]

dt = np.abs(tow - nav_satelity[:, 17])
indeks_najmniejszej_roznicy = np.argmin(dt)

nav_wybrane = nav_satelity[indeks_najmniejszej_roznicy, :]