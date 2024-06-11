import matplotlib.pyplot as plt
import numpy as np
from fun import Hirvonen, Neu

CBKA = np.genfromtxt('projekt3\wyniki\CBKA.pos', comments='%') # wektor krótki
MIMA = np.genfromtxt('projekt3\wyniki\MIMA.pos', comments='%') # wektor średni
WAT = np.genfromtxt('projekt3\wyniki\WAT.pos', comments='%') # wektor bardzo długi

xyzref = [3655333.847, 1403901.067, 5018038.047]
blhref = Hirvonen(xyzref[0], xyzref[1], xyzref[2])
neu = Neu(blhref[0], blhref[1])


