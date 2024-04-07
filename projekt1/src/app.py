import streamlit as st
from funcs import flh2xyz, get_gps_time, Rneu
from alm_module import get_alm_data_str, create_prn_alm2
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

def add_sat_data(t, sat, az, el, visible, x, y, z):
    prn = create_prn_alm2(sat)
    system = check_gnss(prn)
    delta = datetime.timedelta(seconds=t)
    t = startdate + delta
    if t not in satelite_data:
        satelite_data[t] = {}
    if system not in satelite_data[t]:
        satelite_data[t][system] = {}
    satelite_data[t][system][prn] = {'azymut': az, 'elewacja': el, 'visible': visible, 'x': x, 'y': y, 'z': z}

def add_dop_data(t, GDOP, PDOP, HDOP, VDOP, TDOP):
    delta = datetime.timedelta(seconds=t)
    t = startdate + delta
    if t not in dop_data:
        dop_data[t] = {'GDOP': GDOP, 'PDOP': PDOP, 'HDOP': HDOP, 'VDOP': VDOP, 'TDOP': TDOP}

def get_sats_of_chosen_systems():
    sats = df['satelita'].unique()
    sats_of_chosen_systems = []
    for sat in sats:
        system = check_gnss(sat)
        if system == 'GPS' and gps_check == False:
            continue
        if system == 'GLONASS' and glonass_check == False:
            continue
        if system == 'Galileo' and galileo_check == False:
            continue
        if system == 'BeiDou' and beidou_check == False:
            continue
        if system == 'QZSS' and qzss_check == False:
            continue
        if system == 'SBAS' and sbas_check == False:
            continue
        sats_of_chosen_systems.append(sat)
    return sats_of_chosen_systems

def get_chosen_systems():
    systems = []
    if gps_check:
        systems.append('GPS')
    if glonass_check:
        systems.append('GLONASS')
    if galileo_check:
        systems.append('Galileo')
    if beidou_check:
        systems.append('BeiDou')
    if qzss_check:
        systems.append('QZSS')
    if sbas_check:
        systems.append('SBAS')
    return systems

def dop_category(gdop):
    if gdop < 1:
        return 'idealne'
    elif gdop < 2:
        return 'znakomite'
    elif gdop < 4:
        return 'dobre'
    elif gdop < 7:
        return 'umiarkowane'
    elif gdop < 9:
        return 'słabe'
    else:
        return 'złe'

def calculate_dops_for_selected_systems():
    dop_data = {}
    chosen_systems = get_chosen_systems()
    if (len(chosen_systems) == 0):
        st.warning('wybierz co najmniej jeden z systemów GNSS')
        return [], [], [], [], [], [], []
    for t in range (sow0, sow0 + period * 60 * 60, interval*60):
        A = []
        As = []
        for sat in nav:
            prn = create_prn_alm2(sat)
            if check_gnss(prn) not in chosen_systems:
                continue
            x, y, z = satpos(t, week0, sat)
            sat_odb = np.array([x, y, z]) - np.array([odbiornik_x, odbiornik_y, odbiornik_z])
            neu = Rneu(FIrad, LAMBDArad).T.dot(sat_odb)
            az = math.degrees(math.atan2(neu[1], neu[0]))
            if az < 0:
                az += 360
            ro = np.sqrt(neu[0]**2 + neu[1]**2 + neu[2]**2)
            el = math.degrees(math.asin(neu[2]/ro))
            if el > maska:
                a = [-(x-odbiornik_x)/ro, -(y-odbiornik_y)/ro, -(z-odbiornik_z)/ro, 1]
                A.append(a)
        As = np.array(A)
        Q = np.linalg.inv(As.T.dot(As))
        Qdiag = np.diag(Q)
        
        try:
            GDOP = math.sqrt(Qdiag[0] + Qdiag[1] + Qdiag[2] + Qdiag[3])
            PDOP = math.sqrt(Qdiag[0] + Qdiag[1] + Qdiag[2])
            TDOP = math.sqrt(Qdiag[3])

            Qneu = Rneu(math.radians(FI), math.radians(LAMBDA)).T.dot(Q[0:3,0:3]).dot(Rneu(FIrad, LAMBDArad))
            Qneudiag = np.diag(Qneu)

            HDOP = math.sqrt(Qneudiag[0] + Qneudiag[1])
            VDOP = math.sqrt(Qneudiag[2])
        except:
            continue
        delta = datetime.timedelta(seconds=t)
        time = startdate + delta - datetime.timedelta(seconds=sow0)
        if time not in dop_data:
            dop_data[time] = {'GDOP': GDOP, 'PDOP': PDOP, 'HDOP': HDOP, 'VDOP': VDOP, 'TDOP': TDOP}
    time = []
    GDOP = []
    PDOP = []
    HDOP = []
    VDOP = []
    TDOP = []
    category = []
    for t, dops in dop_data.items():
        time.append(t)
        GDOP.append(dops['GDOP'])
        PDOP.append(dops['PDOP'])
        HDOP.append(dops['HDOP'])
        VDOP.append(dops['VDOP'])
        TDOP.append(dops['TDOP'])
        category.append(dop_category(dops['GDOP']))
    return time, GDOP, PDOP, HDOP, VDOP, TDOP, category

st.title('sns - projekt 1 :satellite:')
st.sidebar.header('parametry wejściowe')
date = st.sidebar.date_input('data początkowa:', format='DD.MM.YYYY', value=datetime.date(2024, 2, 29))
time = st.sidebar.time_input('godzina początkowa:', step=3600, value=datetime.time(0, 0))
latitude = st.sidebar.number_input('szerokość geograficzna odbiornika [°]:', min_value=-90.0, max_value=90.0, format="%.6f", step=1.0, value=52.0)
longitude = st.sidebar.number_input('długość geograficzna odbiornika [°]:', min_value=-180.0, max_value=180.0, format="%.6f", step=1.0, value=21.0)
height = st.sidebar.number_input('wysokość odbiornika [m]:', min_value=0.0, format="%.2f", step=1.0, value=100.0)
period = st.sidebar.number_input('długość obserwacji [h]:', min_value=1, max_value=24, step=6, value=24)
interval = st.sidebar.number_input('interwał [min]:', min_value=1, max_value=60, step=1, value=10)
mask = st.sidebar.number_input('maska elewacji [°]:', min_value=0.0, format="%.2f", step=1.0, value=10.0)
gps_check = st.sidebar.checkbox('GPS', value=True)
glonass_check = st.sidebar.checkbox('GLONASS', value=True, key='GLONASS')
galileo_check = st.sidebar.checkbox('Galileo', value=True, key='Galileo')
beidou_check = st.sidebar.checkbox('BeiDou', value=True, key='BeiDou')
qzss_check = st.sidebar.checkbox('QZSS', value=True, key='QZSS')
sbas_check = st.sidebar.checkbox('SBAS', value=True, key='SBAS')
choose = st.selectbox('wybierz wykres:', ['wykres liniowy elewacji', 'wykres liczby widocznych satelitów', 'wykresy DOPów', 'skyplot animacja', 'wykres 3d ruchu satelity'])

MI = 3.986005e14
OMEGA_E = 7.2921151467e-5

# file = r'projekt1/Almanac2024053.alm'
file = r'projekt1\Almanac2024053.alm'

FI = latitude
LAMBDA = longitude
H = height
maska = mask
startdate = datetime.datetime.combine(date, time)
interval = interval
period = period

FIrad = math.radians(FI)
LAMBDArad = math.radians(LAMBDA)
odbiornik_x, odbiornik_y, odbiornik_z = flh2xyz(FIrad, LAMBDArad, H)

nav, prn = get_alm_data_str(file)
satelity = nav[:,0]<500
nav = nav[satelity,:] 
prn = np.array(prn)[satelity]
wiersz_nav = nav[0,:]

week0, sow0 = get_gps_time(startdate.year, startdate.month, startdate.day, startdate.hour, startdate.minute, startdate.second)
week0 = int(week0)
sow0 = int(sow0)

A = []
satelite_data = {}
dop_data = {}
lista_danych = []

for t in range (sow0, sow0 + period * 60 * 60, interval*60):
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
        
        add_sat_data(t-sow0, sat, az, el, visible, x, y, z)
    
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

hover_style = dict(
    font_size=16,
    font_family='Candara'
)

if choose == 'wykres liniowy elewacji':
    elev = go.Figure()
    chosen_systems = get_chosen_systems()
    if chosen_systems == []:
        st.warning('wybierz co najmniej jeden z systemów GNSS!')
        st.stop()
    for system in df['system'].unique():
        if system not in chosen_systems:
            continue
        system_data = df[df['system'] == system]
        for sat in system_data['satelita'].unique():
            sat_data = system_data[system_data['satelita'] == sat]
            elev.add_trace(go.Scatter(
                x=sat_data['czas'],
                y=sat_data['elewacja'],
                mode='lines',
                name=f'{sat}',
                hovertemplate='czas: %{x} <br>elewacja: %{y}°',
                hoverlabel=hover_style,
            ))
    elev.update_layout(
        scene=dict(
            xaxis=dict(title='czas', range=[0, 24]),
            yaxis=dict(title='elewacja'),
        ),
        showlegend=True,
        yaxis_range=[maska, 90],
        yaxis_title='elewacja [°]',
        title='elewacja satelitów w zależności od czasu'
    )
    st.plotly_chart(elev)
elif choose == 'wykres liczby widocznych satelitów':
    visible_sat = df[df['visible'] == True]
    visible_sat = visible_sat.groupby(['czas', 'system']).agg(liczba_widocznych=('visible', 'sum'))
    visible_sat = visible_sat.reset_index()
    chosen_systems = get_chosen_systems()
    if chosen_systems == []:
        st.warning('wybierz co najmniej jeden z systemów GNSS!')
        st.stop()
    if gps_check == False:
        visible_sat = visible_sat[visible_sat['system'] != 'GPS']
    if glonass_check == False:
        visible_sat = visible_sat[visible_sat['system'] != 'GLONASS']
    if galileo_check == False:
        visible_sat = visible_sat[visible_sat['system'] != 'Galileo']
    if beidou_check == False:
        visible_sat = visible_sat[visible_sat['system'] != 'BeiDou']
    if qzss_check == False:
        visible_sat = visible_sat[visible_sat['system'] != 'QZSS']
    if sbas_check == False:
        visible_sat = visible_sat[visible_sat['system'] != 'SBAS']
    vis = px.bar(visible_sat, x='czas', y='liczba_widocznych', color="system", title='liczba widocznych satelitów w zależności od czasu', labels={'czas': 'czas', 'liczba_widocznych': 'liczba widocznych satelitów'})
    vis.update_layout(
        hoverlabel=hover_style
    )
    st.plotly_chart(vis)
elif choose == 'wykresy DOPów':
    time, GDOP, PDOP, HDOP, VDOP, TDOP, category = calculate_dops_for_selected_systems()
    if len(time) == 0:
        st.stop()
    dop = go.Figure()
    dop.add_trace(go.Scatter(x=time, y=GDOP, mode='lines', name='GDOP'))
    dop.add_trace(go.Scatter(x=time, y=PDOP, mode='lines', name='PDOP'))
    dop.add_trace(go.Scatter(x=time, y=HDOP, mode='lines', name='HDOP'))
    dop.add_trace(go.Scatter(x=time, y=VDOP, mode='lines', name='VDOP'))
    dop.add_trace(go.Scatter(x=time, y=TDOP, mode='lines', name='TDOP'))
    dop.update_layout(title='wykres liniowy DOPów w zależności od czasu', yaxis_title='DOP', hoverlabel=hover_style)
    st.plotly_chart(dop)

    cat = px.pie(names=['idealne', 'znakomite', 'dobre', 'umiarkowane', 'słabe', 'złe'], values=[category.count('idealne'), category.count('znakomite'), category.count('dobre'), category.count('umiarkowane'), category.count('słabe'), category.count('złe')], title='warunki pomiaru', labels={'names': 'kategoria', 'values': 'liczba'}, color_discrete_sequence=px.colors.sequential.Aggrnyl)
    st.plotly_chart(cat)

    st.markdown('**przyjęte przedziały:**')
    st.caption('idealne: GDOP < 1')
    st.caption('znakomite: 1 <= GDOP < 2')
    st.caption('dobre: 2 <= GDOP < 4')
    st.caption('umiarkowane: 4 <= GDOP < 7')
    st.caption('słabe: 7 <= GDOP < 9')
    st.caption('złe: GDOP >= 9')
elif choose == 'skyplot animacja':
    selected_systems = get_chosen_systems()
    if selected_systems == []:
        st.warning('wybierz co najmniej jeden z systemów GNSS!')
        st.stop()
    filtered_df = df[df['system'].isin(selected_systems)]
    skyplot = px.scatter_polar(
        filtered_df, r='elewacja', theta='azymut',
        color='system', animation_frame='czas', range_r=[maska, 90],
        hover_data={'czas': False, 'system': True, 'satelita': False},
        hover_name='satelita'
    )
    st.plotly_chart(skyplot)
elif choose == 'wykres 3d ruchu satelity':
    r = 6378137
    u, v = np.mgrid[0:(2 * np.pi+0.1):0.2, -np.pi/2:np.pi/2:0.2]
    x = r*np.cos(v) * np.cos(u)
    y = r*np.cos(v) * np.sin(u)
    z = r*np.sin(v)
    plot3d = go.Figure()
    plot3d.add_trace(go.Surface(x=x, y=y, z=z, colorscale='blues'))
    plot3d.update_layout(scene=dict(aspectratio=dict(x=1, y=1, z=1),
                                aspectmode='manual'))
    
    sats = get_sats_of_chosen_systems()
    if len(sats) == 0:
        st.warning('wybierz co najmniej jeden z systemów GNSS!')
        st.stop()
    satelite_choice = st.selectbox('wybierz satelitę:', sats)
    sat_data = df[df['satelita'] == satelite_choice]
    plot3d.add_trace(go.Scatter3d(x=sat_data['x'], y=sat_data['y'], z=sat_data['z'], mode='lines', name=f'{satelite_choice}'))
    plot3d.update_layout(title=f'wykres 3d ruchu satelity {satelite_choice} w układzie ECEF', hoverlabel=hover_style)
    st.plotly_chart(plot3d)
