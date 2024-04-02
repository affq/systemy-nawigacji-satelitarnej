import streamlit as st
import projekt1 as p1
from streamlit_option_menu import option_menu

def on_change(key):
    selection = st.session_state[key]


selected2 = option_menu("systemy nawigacji satelitarnej - projekt 1", ["wykres elewacji", "widoczne satelity", "wykres dop-ów", 'skyplot'], 
icons=['align-bottom', 'radar', 'graph-up', 'pie-chart'], 
default_index=0, orientation="horizontal", on_change=on_change,
key="mainmenu")

st.sidebar.header('Parametry')

date = st.sidebar.date_input('Wybierz datę:', value=p1.startdate, format='DD.MM.YYYY')
time = st.sidebar.time_input('Wybierz godzinę:', step=3600, value=p1.startdate.time())

latitude = st.sidebar.number_input('Szerokość geograficzna:', value=float(p1.FI), min_value=-90.0, max_value=90.0, format="%.6f")
longitude = st.sidebar.number_input('Długość geograficzna:', value=float(p1.LAMBDA), min_value=-180.0, max_value=180.0, format="%.6f")
height = st.sidebar.number_input('Wysokość:', value=float(p1.H), min_value=0.0, format="%.2f")
