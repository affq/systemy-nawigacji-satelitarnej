import numpy as np
import math

a = 6378137 # wielka półoś elipsoidy GRS80 w metrach
e2 = 0.00669438002290 #kwadrat pierwszego mimośrodu dla elipsoidy GRS80

def julday(y,m,d,h=0):
    if m <= 2:
        y = y - 1
        m = m + 12
    jd = math.floor(365.25*(y+4716))+math.floor(30.6001*(m+1))+d+h/24-1537.5
    return jd

def get_gps_time(y,m,d,h=0,mnt=0,s=0):
    days = julday(y,m,d) - julday(1980,1,6)
    week = days//7
    day = days%7
    sow = day * 86400 + h * 3600 + mnt * 60 + s
    return week, sow

def flh2xyz(phi, lamb, h):
    '''
    funkcja przeliczająca współrzędne geograficzne na współrzędne kartezjańskie

    phi, lamb  w radianach
    '''
    N = a/np.sqrt(1-e2*np.sin(phi)**2)
    x = (N + h)*np.cos(phi)*np.cos(lamb)
    y = (N + h)*np.cos(phi)*np.sin(lamb)
    z = (N*(1-e2)+h)*np.sin(phi)
    return [x, y, z]

def Rneu(phi, lamb):
    '''
    funkcja zwracająca macierz rotacji z układu ECEF do układu NEU

    phi, lamb  w radianach
    '''
    R = np.array([[-np.sin(phi)*np.cos(lamb), -np.sin(lamb), np.cos(phi)*np.cos(lamb)],
                    [-np.sin(phi)*np.sin(lamb), np.cos(lamb), np.cos(phi)*np.sin(lamb)],
                    [np.cos(phi), 0, np.sin(phi)]])
    return R