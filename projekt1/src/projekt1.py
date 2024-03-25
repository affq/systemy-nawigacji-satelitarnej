import numpy as np
from alm_module import get_alm_data_str, get_prn_number2, create_prn_alm2
from funcs import get_gps_time, flh2xyz, Rneu
import math
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import random

FI = 52 #szerokość geograficzna odbiornika
LAMBDA = 21 #długość geograficzna odbiornika
H = 100 #wysokość odbiornika w metrach

FIrad = math.radians(FI)
LAMBDArad = math.radians(LAMBDA)

year, m, d = 2024, 2, 29 # czas na który chcemy wyznaczyć pozycję satelitów
h, mnt, s = 12, 0, 0

# +18s bo sekundy przestępne

MI = 3.986005e14
OMEGA_E = 7.2921151467e-5

file = 'projekt1\Almanac2024053.alm'

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

odbiornik_x, odbiornik_y, odbiornik_z = flh2xyz(FIrad, LAMBDArad, H) # współrzędne odbiornika w układzie xyz

print(f"Wyniki dla epoki y={year} m={m} d={d} h={h} mnt={mnt} s={s}")
print(f"Dane z pliku almanach: {file} \n")

print(f'Współrzędne odbiornika phi: {FI}, lambda: {LAMBDA} [deg]; h: {h}m')
print(f'Współrzędne XYZ odbiornika: {odbiornik_x}, {odbiornik_y}, {odbiornik_z}')
print(f'Maska elewacji: {maska} [deg] \n')

# obliczenie co 10 minut pozycji satelitow
h0, mnt0, s0 = 0, 0, 0
week0, sow0 = get_gps_time(year,m,d,h0,mnt0,s0)

#week0 and sow0 to int
week0 = int(week0)
sow0 = int(sow0)

satelite_data = {}
dop_data = {}

def add_sat_data(t, sat, az, el, visible):
    prn = create_prn_alm2(sat)
    if t not in satelite_data:
        satelite_data[t] = {}
    if prn not in satelite_data[t]:
        satelite_data[t][prn] = {'azimuth': az, 'elevation': el, 'visible': visible}

def add_dop_data(t, GDOP, PDOP, HDOP, VDOP, TDOP):
    if t not in dop_data:
        dop_data[t] = {'GDOP': GDOP, 'PDOP': PDOP, 'HDOP': HDOP, 'VDOP': VDOP, 'TDOP': TDOP}

for t in range (sow0, sow0 + 24 * 60 * 60, 600):
    for sat in nav:
        x, y, z = satpos(t, week0, sat) # współrzędne satelity w układzie xyz
        sat_odb = np.array([x, y, z]) - np.array([odbiornik_x, odbiornik_y, odbiornik_z]) # wektor satelita-odbiornik xyz
        neu = Rneu(FIrad, LAMBDArad).T.dot(sat_odb) # wektor satelita-odbiornik neu
        az = math.degrees(math.atan2(neu[1], neu[0])) # azymut [deg]
        print(f"SATELITA: {int(sat[0])}")
        if az < 0:
            az += 360
        print(f"Azymut: {round(az, 3)} [deg]")
        ro = np.sqrt(neu[0]**2 + neu[1]**2 + neu[2]**2) # odległość [m]
        el = math.degrees(math.asin(neu[2]/ro)) # elewacja [deg]
        print(f"Elewacja: {round(el, 3)} [deg]")

        if el > maska:
            a = [-(x-odbiornik_x)/ro, -(y-odbiornik_y)/ro, -(z-odbiornik_z)/ro, 1]
            A.append(a)
            visible = True
        else:
            visible = False
        
        add_sat_data((t-sow0)/3600, sat, az, el, visible)
    
    
    As = np.array(A)
    Q = np.linalg.inv(As.T.dot(As))
    Qdiag = np.diag(Q)

    GDOP = math.sqrt(Qdiag[0] + Qdiag[1] + Qdiag[2] + Qdiag[3])
    PDOP = math.sqrt(Qdiag[0] + Qdiag[1] + Qdiag[2])
    TDOP = math.sqrt(Qdiag[3])

    Qneu = Rneu(math.radians(FI), math.radians(LAMBDA)).T.dot(Q[0:3,0:3]).dot(Rneu(FIrad, LAMBDArad))
    Qneudiag = np.diag(Qneu)

    HDOP = math.sqrt(Qneudiag[0] + Qneudiag[1])
    VDOP = math.sqrt(Qneudiag[2])

    add_dop_data((t-sow0)/3600, GDOP, PDOP, HDOP, VDOP, TDOP)

# wykres liniowy elewacji - godziny x elewacje y;
# dobrym pomysłem usuwać łuki poniżej maski obserwacji
# line chart

df = pd.DataFrame.from_dict(satelite_data)
df = df.stack().apply(pd.Series).reset_index()
df.columns = ['satelita', 'czas', 'azymut', 'elewacja', 'widoczność']
print(df)

fig = go.Figure()
for sat in df['satelita'].unique():
    sat_data = df[df['satelita'] == sat]
    fig.add_trace(go.Scatter(
        x=sat_data['czas'],
        y=sat_data['elewacja'],
        mode='lines',
        name=f'Satelita {sat}'
    ))

fig.update_layout(
    scene=dict(
        xaxis=dict(title='czas', range=[0, 24]),
        yaxis=dict(title='elewacja', range=[0, 90]),
    ),
    title='Pozycje satelitów',
    showlegend=True
)

fig.show()
exit()

# wykres z liczbą widocznych satelitów (elewacja większa od maski) w zależności od czasu
# może wykres słupkowy lepiej
# bar chart
fig = go.Figure()
fig = px.bar(x=list(satelite_data.keys()), y=[sum([1 for sat in satelite_data[time].values() if sat['visible']]) for time in dop_data.keys()], labels={'x': 'Czas [h]', 'y': 'Liczba widocznych satelitów'})
fig.show()

# liniowy w zaleznosci od wartosci dopow
# line chart


# kolowy w zaleznosci od wartosci procentowej dop w ciagu doby

# wykres skyplot na wybraną godzinę
# łuków nie trzeba ale położenie tak
    
#mapka z ruchem satelitów 