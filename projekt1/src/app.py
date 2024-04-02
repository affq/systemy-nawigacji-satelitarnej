import streamlit as st
import projekt1 as p1
from funcs import flh2xyz, get_gps_time, Rneu
from alm_module import get_alm_data_str, create_prn_alm2
from streamlit_option_menu import option_menu
import plotly.express as px
import plotly.graph_objects as go
import datetime
import math
import numpy as np
import pandas as pd

def check_gnss(prn):
    system = ''
    if prn[0] == 'G':
        system = 'GPS'
    elif prn[0] == 'R':
        system = 'GLONASS'
    elif prn[0] == 'E':
        system = 'Galileo'
    elif prn[0] == 'C':
        system = 'BeiDou'
    elif prn[0] == 'Q':
        system = 'QZSS'
    elif prn[0] == 'S':
        system = 'SBAS'
    else:
        system = 'Unknown'
    return system

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

def add_sat_data(t, sat, az, el, visible):
    prn = create_prn_alm2(sat)
    system = check_gnss(prn)
    delta = datetime.timedelta(seconds=t)
    t = startdate + delta
    if t not in satelite_data:
        satelite_data[t] = {}
    if system not in satelite_data[t]:
        satelite_data[t][system] = {}
    satelite_data[t][system][prn] = {'azymut': az, 'elewacja': el, 'visible': visible}

def add_dop_data(t, GDOP, PDOP, HDOP, VDOP, TDOP):
    delta = datetime.timedelta(seconds=t)
    t = startdate + delta
    if t not in dop_data:
        dop_data[t] = {'GDOP': GDOP, 'PDOP': PDOP, 'HDOP': HDOP, 'VDOP': VDOP, 'TDOP': TDOP}

MI = 3.986005e14
OMEGA_E = 7.2921151467e-5

file = r'C:\Users\adria\Desktop\systemy-nawigacji-satelitarnej\projekt1\Almanac2024053.alm'

FI = 52 #szerokość geograficzna odbiornika
LAMBDA = 21 #długość geograficzna odbiornika
H = 100 #wysokość odbiornika w metrach
maska = 10 # maska elewacji w stopniach
year, m, d = 2024, 2, 29
h0, mnt0, s0 = 0, 0, 0
startdate = datetime.datetime(year, m, d, h0, mnt0, s0)
# +18s bo sekundy przestępne

FIrad = math.radians(FI)
LAMBDArad = math.radians(LAMBDA)
odbiornik_x, odbiornik_y, odbiornik_z = flh2xyz(FIrad, LAMBDArad, H) # współrzędne odbiornika w układzie xyz

A = []

nav, prn = get_alm_data_str(file)
satelity = nav[:,0]<500
nav = nav[satelity,:] 
prn = np.array(prn)[satelity]
wiersz_nav = nav[0,:]

week0, sow0 = get_gps_time(year,m,d,h0,mnt0,s0)
week0 = int(week0)
sow0 = int(sow0)

satelite_data = {}
dop_data = {}

for t in range (sow0, sow0 + 24 * 60 * 60, 600):
    A = []
    As = []
    for sat in nav:
        x, y, z = satpos(t, week0, sat) # współrzędne satelity w układzie xyz
        sat_odb = np.array([x, y, z]) - np.array([odbiornik_x, odbiornik_y, odbiornik_z]) # wektor satelita-odbiornik xyz
        neu = Rneu(FIrad, LAMBDArad).T.dot(sat_odb) # wektor satelita-odbiornik neu
        az = math.degrees(math.atan2(neu[1], neu[0])) # azymut [deg]
        if az < 0:
            az += 360
        ro = np.sqrt(neu[0]**2 + neu[1]**2 + neu[2]**2) # odległość [m]
        el = math.degrees(math.asin(neu[2]/ro)) # elewacja [deg]

        if el > maska:
            a = [-(x-odbiornik_x)/ro, -(y-odbiornik_y)/ro, -(z-odbiornik_z)/ro, 1]
            A.append(a)
            visible = True
        else:
            visible = False
        
        add_sat_data(t-sow0, sat, az, el, visible)
    
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

    add_dop_data(t-sow0, GDOP, PDOP, HDOP, VDOP, TDOP)


lista_danych = []



for czas, systemy in satelite_data.items():
    for system, satelity in systemy.items():
        for satelita, dane in satelity.items():
            dane = {
                "czas": czas,
                "system": system,
                "satelita": satelita,
                **dane
            }
            lista_danych.append(dane)

df = pd.DataFrame(lista_danych)

elev = go.Figure()
for sat in df['satelita'].unique():
    sat_data = df[df['satelita'] == sat]
    elev.add_trace(go.Scatter(
        x=sat_data['czas'],
        y=sat_data['elewacja'],
        mode='lines',
        name=f'Satelita {sat}'
    ))

elev.update_layout(
    scene=dict(
        xaxis=dict(title='czas', range=[0, 24]),
        yaxis=dict(title='elewacja'),
    ),
    title='Pozycje satelitów',
    showlegend=True,
    yaxis_range=[maska, 90]
)

visible_sat = df[df['visible'] == True]
visible_sat = visible_sat.groupby(['czas', 'system']).agg(liczba_widocznych=('visible', 'sum'))
visible_sat = visible_sat.reset_index()
vis = px.bar(visible_sat, x='czas', y='liczba_widocznych', color="system", title='Liczba widocznych satelitów w zależności od czasu') 


dops_df = pd.DataFrame(dop_data).T.reset_index()
dops_df = dops_df.rename(columns={'index': 'czas'})
lines = []
for dop_type in ['GDOP', 'PDOP', 'HDOP', 'VDOP', 'TDOP']:
    lines.append(go.Scatter(x=dops_df['czas'], y=dops_df[dop_type], mode='lines', name=dop_type))
dop = go.Figure(lines)
dop.update_layout(title='Wykres DOP', xaxis_title='Czas', yaxis_title='Wartość')


skyplot = px.scatter_polar(
    df, r='elewacja', theta='azymut',
    color='system', animation_frame='czas', range_r=[maska, 90]
)

def show_chart(key):
    if key == 'wykres elewacji':
        st.plotly_chart(elev)
    elif key == 'widoczne satelity':
        st.plotly_chart(vis)
    elif key == 'wykres dop-ów':
        st.plotly_chart(dop)
    elif key == 'skyplot':
        st.plotly_chart(skyplot)

st.title('Systemy nawigacji satelitarnej - projekt 1')

choose = st.selectbox('Wybierz wykres:', ['wykres elewacji', 'widoczne satelity', 'wykres dop-ów', 'skyplot'])
print(choose)
show_chart(choose)

st.sidebar.header('Parametry wejściowe')

date = st.sidebar.date_input('Data :', value=p1.startdate, format='DD.MM.YYYY')
time = st.sidebar.time_input('Godzina:', step=3600, value=p1.startdate.time())

latitude = st.sidebar.number_input('Szerokość geograficzna [°]:', value=float(p1.FI), min_value=-90.0, max_value=90.0, format="%.6f", step=1.0)
longitude = st.sidebar.number_input('Długość geograficzna [°]:', value=float(p1.LAMBDA), min_value=-180.0, max_value=180.0, format="%.6f", step=1.0)
height = st.sidebar.number_input('Wysokość [m]:', value=float(p1.H), min_value=0.0, format="%.2f", step=1.0)
mask = st.sidebar.number_input('Maska elewacji [°]:', value=float(p1.maska), min_value=0.0, format="%.2f", step=1.0)

submit = st.sidebar.button('Oblicz', use_container_width=True, on_click=lambda: st.experimental_rerun())
