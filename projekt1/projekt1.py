import wget
import alm_module as alm

file = 'ftp://ftp.trimble.com/pub/eph/Almanac.alm'
wget.download(file, 'Almanac.alm')