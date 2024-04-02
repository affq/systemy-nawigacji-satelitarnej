import streamlit as st
import projekt1 as p1
from projekt1 import FI, LAMBDA, H, maska
from streamlit_option_menu import option_menu

def on_change(key):
    selection = st.session_state[key]

def calculate():
    return 0

selected2 = option_menu("systemy nawigacji satelitarnej - projekt 1", ["wykres elewacji", "widoczne satelity", "wykres dop-ów", 'skyplot'], 
icons=['align-bottom', 'radar', 'graph-up', 'pie-chart'], 
default_index=0, orientation="horizontal", on_change=on_change,
key="mainmenu")

st.sidebar.header('Parametry wejściowe')

date = st.sidebar.date_input('Data :', value=p1.startdate, format='DD.MM.YYYY')
time = st.sidebar.time_input('Godzina:', step=3600, value=p1.startdate.time())

latitude = st.sidebar.number_input('Szerokość geograficzna [°]:', value=float(p1.FI), min_value=-90.0, max_value=90.0, format="%.6f", step=1.0)
longitude = st.sidebar.number_input('Długość geograficzna [°]:', value=float(p1.LAMBDA), min_value=-180.0, max_value=180.0, format="%.6f", step=1.0)
height = st.sidebar.number_input('Wysokość [m]:', value=float(p1.H), min_value=0.0, format="%.2f", step=1.0)
mask = st.sidebar.number_input('Maska elewacji [°]:', value=float(p1.maska), min_value=0.0, format="%.2f", step=1.0)

submit = st.sidebar.button('Oblicz', use_container_width=True, on_click=calculate)
