import numpy as np
from readrnx_studenci import readrnxnav, readrnxobs, date2tow
from funcs import *

OMEGA = 7.2921151467e-5 #[rad/s]
C = 299792458.0 #[m/s]

# cieżka do pliku nawigacyjnego
nav_file = 'rinex/BRDC00WRD_R_20240650000_01D_GN.rnx'
# cieżka do pliku obserwacyjnego
obs_file = 'rinex/JOZ200POL_R_20240650000_01D_30S_MO.rnx'

# zdefiniowanie czasu obserwacji: daty początkowej i końcowej
# dla pierwszej epoki z pliku będzie to:
time_start = [2024, 3, 5, 0, 0, 0]  
time_end = [2024, 3, 5, 23, 59, 59] 

# odczytanie danych z pliku obserwacyjnego
obs, iobs = readrnxobs(obs_file, time_start, time_end, 'G')
# odczytanie danych z pliku nawigacyjnego:
nav, inav = readrnxnav(nav_file)

# print(obs, iobs)
# print(nav, inav)

"""
zdefiniowanie współrzędnych przybliżonych odbiornika - mogą to być współrzędne z nagłówka 
pliku obserwacyjnego, skopiowane "z palca" lub pobierane automatycznie z treci nagłówka pliku Rinex
"""
xr0 = [3660000.,  1400000.,  5000000.]

"""
Wprowadzenie ustawień, takich jak maska obserwacji, czy typ poprawki troposferycznej
"""
el_mask = 10 # elevation mask/cut off in degrees
 
"""
Przeliczenie daty początkowej i końcowej do sekund tygodnia GPS - niezbędne w celu
poprawnej definicji pętli związanej z czasem obserwacji w ciągu całej doby
"""
week, tow = date2tow(time_start)[0:2]
week_end, tow_end = date2tow(time_end)[0:2]

"""
Otwieramy dużą pętlę
for t in range(tow, tow_end+1, dt): gdzie dt równe 30
"""
t = 213300

index_t = iobs[:,2]==t
Pobs = obs[index_t,0]
satelity = iobs[index_t,0]

tau = 0.07
dtr = 0

"""
Wewnątrz tej pętli, zajmujemy się obserwacjami wyłącznie dla jednej epoki (epoka t), zatem:
    1. Wybieramy obserwacje dla danej epoki, na podstawie tablicy iobs oraz naszej epoki t
    czyli, wybieramy te obserwacje z tablicy obs, dla których w tablicy iobs ostatnia kolumna 
    jest równa t - przypisujemy do zmiennej np. Pobs
    2. wybieramy satelity, obserwowane w danej epoce, na podstawie tablicy iobs - na podstawie 
    naszego t - przypisujemy do zmiennej np. sats
    3. Definiujemy wartości przybliżone błąd zegara odbiornika
    dtr = 0 oraz czasu propagacji sygnału tau = 0.07
    4. Najprawdopodobniej przyda się definicja pustych wektorów, np. zawierających 
    odległosci geometryczne (wartoci przybliżone na podstawie tau)

        
    Przechodzimy do iteracyjnego obliczenia współrzędnych odbiornika - w pierwszych testach naszego programu, zróbmy obliczenia nieiteracyjnie, 
    ale pamiętajmy o tym, że będzie trzeba przygotować kod do działania w pętli:
        
        Po weryfikacji działania programu, można zamienić pętlę for na pętle while, dopisując
        warunek zbieżnoci kolejnych współrzędnych - skróci nam to czas obliczeń, ponieważ 
        najczęściej wystarcza mniej iteracji niż 5
        """
for i in range(5):
    for k, sat in enumerate(satelity):
        ts = t - tau + dtr
        alfa = tau * OMEGA
        R = np.array([[np.cos(alfa), np.sin(alfa), 0], [-np.sin(alfa), np.cos(alfa), 0], [0, 0, 1]])
        xs_rot = R@xr0
        rho = np.sqrt((xs_rot[0] - xr0[0])**2 + (xs_rot[1] - xr0[1])**2 + (xs_rot[2] - xr0[2])**2)
        """
                    Obliczamy elewację i azymut
                    Macierz Rneu definiujemy na podstawie x0, przeliczonego do współrzędnych
                    phi, lambda, algorytmem Hirvonena
        """
        B, L, H = hirvonen(xr0[0], xr0[1], xr0[2])
        RneuMatrix = Rneu(B, L)
        """
                    Odrzucamy satelity znajdujące się poniżej maski
                    
                        Obliczamy poprawki atmosferyczne - dopiero wówczas, kiedy działać będzie nam program bez uwzględniania poprawek:
                            trop oraz iono
                    
                    Wyznaczamy pseudoodległosć przybliżoną (obliczoną), jako:
                        Pcalc = rho - cdts + dtr + trop + iono
                        
                    Wyznaczamy kolejne elementy wektora wyrazów wolnych y, jako:
                        y = Pobs - Pcalc
                        
                    Budujemy kolejne wiersze macierzy A:
                
                Kończymy pętle dla kolejnych satelitów
                
                1. Łączymy ze sobą elementy wektora wyrazów wolych w jeden wektor
                2. Łączymy ze sobą kolejnę wiersze macierz współczynników A
                3. Rozwiązujemy układ równań, metodą najmniejszych kwadratów
                
                               
                Aktualizujemy wartosci przybliżone o odpowiednie elementy wektora x
                xr[0] = x0[0] + x[0]
                xr[1] = x0[1] + x[1]
                xr[2] = x0[2] + x[2]
                dtr = dtr + x[3]/c 
                
                Tak obliczone wartoci xr oraz dtr stanowią wartoci wejsciowe do kolejnej iteracji, itd 
                do skończenia piątej iteracji lub spełnienia warunku zbieżnoci współrzędncyh
            
            
            Po skończeniu 5. iteracji, zbieramy obliczone współrzędne xr - warto zebrać również
            liczby obserwowanych satelitów, obliczone wartoci współczynników DOP (przynajmniej PDOP)
            
"""








