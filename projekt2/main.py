import numpy as np
from readrnx_studenci import readrnxnav, readrnxobs, date2tow
from funcs import *
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox

OMEGA = 7.2921151467e-5
C = 299792458.0

nav_file = 'rinex/BRDC00WRD_R_20240650000_01D_GN.rnx'
obs_file = 'rinex/JOZ200POL_R_20240650000_01D_30S_MO.rnx'

time_start = [2024, 3, 5, 0, 0, 0]  
time_end = [2024, 3, 5, 23, 59, 59]

obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
nav, inav = readrnxnav(nav_file)
zdrowe = nav[:,30] == 0
nav = nav[zdrowe,:]
inav = inav[zdrowe]

xr0 = xr_reference = [3664880.9100,  1409190.3850,  5009618.2850]

el_mask = 10
 
week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]

satelity = []

tau = 0.07
dtr = 0
ro = 0
dt = 30

XR5 = []
DIFFS = []

def calc(is_tropo=False, is_iono=False):
    global xr0, dtr, satelity, tau, ro, dt, XR5, DIFFS, root, response
    root.destroy()
    if is_tropo:
        response = messagebox.askyesno("Poprawka troposferyczna", "Czy chcesz użyć modelu Hopfielda do obliczenia poprawki troposferycznej? \nJeżeli nie, zostanie użyty model Saastamoinena.")
    for t in range (tow, tow_end+1, 30):
        index_t = iobs[:,2]==t
        Pobs = obs[index_t,0]
        satelity = iobs[index_t,0]
        taus = []
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
                    tropo_hop = Hopfield(H, np.rad2deg(el), i)
                    tropo_sas = Saastamoinen(H, np.rad2deg(el), i)
                    iono = Klobuchar(np.rad2deg(B), np.rad2deg(L), np.rad2deg(el), np.rad2deg(az), ts)           
                    Pcalc = ro + c*dtr - c*dt0s 
                    if is_tropo:
                        if response:
                            Pcalc += tropo_hop
                        else:
                            Pcalc += tropo_sas
                    if is_iono:
                        Pcalc += iono
                    y = Pobs[nr] - Pcalc
                    Y.append(y)
                    a = [-(xs_rot[0] - xr0[0])/ro, -(xs_rot[1] - xr0[1])/ro, -(xs_rot[2] - xr0[2])/ro, 1]
                    A.append(a)

            A = np.array(A)
            Y = np.array(Y)
            x = np.linalg.pinv(A.T@A)@A.T@Y
            xr0 = [xr0[0] + x[0], xr0[1] + x[1], xr0[2] + x[2]]
            dtr = dtr + x[3]/C
            if i == 4:
                XR5.append([t, xr0[0], xr0[1], xr0[2]])
                diff = [xr0[0] - xr_reference[0], xr0[1] - xr_reference[1], xr0[2] - xr_reference[2]]
                DIFFS.append([t, diff[0], diff[1], diff[2]])

    DIFFS = np.array(DIFFS)
    fig, axes = plt.subplots(3, 1, figsize=(10, 10))
    fig.suptitle(f"Różnice między współrzędnymi obliczonymi i referencyjnymi \n tropo: {is_tropo}, iono: {is_iono}")
    std_x = np.std(DIFFS[:,1])
    rms_x = np.sqrt(np.mean(DIFFS[:,1]**2))
    axes[0].plot((DIFFS[:,0] - tow)/3600, DIFFS[:,1], 'r')
    axes[0].text(0.95, 0.95, f"σ = {std_x:.2f}m \n RMS = {rms_x:.2f}m", transform=axes[0].transAxes, ha='right', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    axes[0].grid()
    axes[0].set_ylabel('dx[m]')
    std_y = np.std(DIFFS[:,2])
    rms_y = np.sqrt(np.mean(DIFFS[:,2]**2))
    axes[1].plot((DIFFS[:,0] - tow)/3600, DIFFS[:,2], 'g')
    axes[1].text(0.95, 0.95, f"σ = {std_y:.2f}m \n RMS = {rms_y:.2f}m", transform=axes[1].transAxes, ha='right', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    axes[1].grid()
    axes[1].set_ylabel('dy[m]')
    std_z = np.std(DIFFS[:,3])
    rms_z = np.sqrt(np.mean(DIFFS[:,3]**2))
    axes[2].plot((DIFFS[:,0] - tow)/3600, DIFFS[:,3], 'b')
    axes[2].text(0.95, 0.95, f"σ = {std_z:.2f}m \n RMS = {rms_z:.2f}m", transform=axes[2].transAxes, ha='right', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))
    axes[2].grid()
    axes[2].set_ylabel('dz[m]')
    plt.show()
    fig.savefig(f"wykresy/tropo_{is_tropo}_iono_{is_iono}.png")

root = tk.Tk()
root.title("snsx2")
root.geometry("300x130")

q_label = tk.Label(root, text="Wybierz, które poprawki uwzględnić w obliczeniach:")
q_label.pack()

tropo_var = tk.BooleanVar()
iono_var = tk.BooleanVar()

tropo_checkbutton = tk.Checkbutton(root, text="troposferyczna", variable=tropo_var)
iono_checkbutton = tk.Checkbutton(root, text="jonosferyczna", variable=iono_var)

tropo_checkbutton.pack()
iono_checkbutton.pack()

calc_button = tk.Button(root, text="Oblicz", command=lambda: calc(tropo_var.get(), iono_var.get()))
calc_button.pack()

root.mainloop()










