import matplotlib.pyplot as plt
import numpy as np
from fun import Hirvonen, Rneu

xyzref = [3655333.847, 1403901.067, 5018038.047]
blhref = Hirvonen(xyzref[0], xyzref[1], xyzref[2])
rneu = Rneu(blhref[0], blhref[1])

interval = 30

CBKA = np.genfromtxt('projekt3\wyniki\CBKA.pos', comments='%') # wektor krótki
dxyz_CBKA = CBKA[:,2:5] - xyzref
dneu_CBKA = []

MIMA = np.genfromtxt('projekt3\wyniki\MIMA.pos', comments='%') # wektor średni
dxyz_MIMA = MIMA[:,2:5] - xyzref
dneu_MIMA = []

WAT = np.genfromtxt('projekt3\wyniki\WAT.pos', comments='%') # wektor bardzo długi
dxyz_WAT = WAT[:,2:5] - xyzref
dneu_WAT = []

for i, xyz in enumerate(dxyz_CBKA):
    neu = rneu.T.dot(xyz)
    dneu_CBKA.append([i*interval/3600, neu[0], neu[1], neu[2]])

for i, xyz in enumerate(dxyz_MIMA):
    neu = rneu.T.dot(xyz)
    dneu_MIMA.append([i*interval/3600, neu[0], neu[1], neu[2]])

for i, xyz in enumerate(dxyz_WAT):
    neu = rneu.T.dot(xyz)
    dneu_WAT.append([i*interval/3600, neu[0], neu[1], neu[2]])


# wykres liniowy błędów współrzędnych w zależności od czasu dla CBKA
plt.figure()
plt.plot([d[0] for d in dneu_CBKA], [d[1] for d in dneu_CBKA], label='dN')
plt.plot([d[0] for d in dneu_CBKA], [d[2] for d in dneu_CBKA], label='dE')
plt.plot([d[0] for d in dneu_CBKA], [d[3] for d in dneu_CBKA], label='dU')
plt.xlabel('Czas [h]')
plt.ylabel('Błąd [m]')
plt.title('Błąd współrzędnych dla CBKA')
plt.legend()
plt.grid()
plt.show()

# wykres liniowy błędów współrzędnych w zależności od czasu dla CBKA - oddzielny dla każdej współrzędnej
fig, axes = plt.subplots(3, 1)
axes[0].plot([d[0] for d in dneu_CBKA], [d[1] for d in dneu_CBKA])
axes[0].set_title('dN')
axes[0].set_xlabel('Czas [h]')
axes[0].set_ylabel('Błąd [m]')
axes[0].grid()

axes[1].plot([d[0] for d in dneu_CBKA], [d[2] for d in dneu_CBKA])
axes[1].set_title('dE')
axes[1].set_xlabel('Czas [h]')
axes[1].set_ylabel('Błąd [m]')
axes[1].grid()

axes[2].plot([d[0] for d in dneu_CBKA], [d[3] for d in dneu_CBKA])
axes[2].set_title('dU')
axes[2].set_xlabel('Czas [h]')
axes[2].set_ylabel('Błąd [m]')
axes[2].grid()

plt.tight_layout()
plt.show()
