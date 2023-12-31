# type: ignore
# flake8: noqa
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: fragment}

import math
AU = 149597870 # km
c = 299792.458 # km/s
ts = AU/c # s
tm = ts/60 # min
print(f"Time taken by light to reach Earth from the Sun = {ts:1.4} s = {tm:1.4} min = 8 min 19 s ")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: fragment}

import math
AU = 149597870 # km
aphD = 152100000 # km
perD = 147100000 # km
AU_avg = (aphD+perD)/(2*AU) # Average distance, gives AU, check: print(f"{AU_avg}")
print(f"Distance Earth-Sun at aphelion (most distant) = {aphD:1.3e} km = {aphD/AU:1.4} au")
print(f"Distance Earth-Sun at perihelion (closest) = {perD:1.3e} km = {perD/AU:1.3} au")
print(f"Average distance Earth-Sun = {AU:1.3e} km = 1 au")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
import numpy as np
import math
AU = 149597870 # km
pc = 648000/np.pi*AU # km
print(f"1parsec = {pc:1.3e} km")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
import math
import numpy as np
def arcsec_to_rad(p):
    rad = p/206265
    return rad

p = 0.7681 # as
d_pc = 1/p # pc

p_rad = arcsec_to_rad(p)
d_au = 1/p_rad # au
d_km = d_au*AU # km
print(f"p = {p} as")
print(f"d = 1/p pc = 1/{p} pc = {d_pc:1.1e} pc")
print(f"d = 1/p pc = 206265/p au = 206265/{p} au = {d_au:1.1e} au")
print(f"d = {d_au:1.1e} au * 149597870 km = {d_km:1.1e} km")
print(f"The distance from Proxima Centauri is = {d_pc:1.1e} pc = {d_au:1.1e} au = {d_km:1.1e} km")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: fragment}
import math
AU = 149597870 # km
ang_diameter_arcm = 32 # angular diameter of Sun is about 32 minutes of arc
ang_diameter = ang_diameter_arcm/60 * math.pi/180 
# ang_diameter = ang_diameter_rad
diameter = AU * ang_diameter
print(f"1 au = {AU:1.3e} km")
print(f"angular diameter of the sun =  {ang_diameter:.3e} radians")
print(f"Diameter of Sun = {diameter:.0f} km")
print(f"Radius of Sun = {diameter/2:.0f} km")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: skip}
import scipy
from scipy.constants import pi, au, parsec

F = 1361.0 # W/m/m
r = au # m
L = F*4*pi*r**2

print(f'Power flux: F = {F} W/m/m')
print(f'r = 1au = {r:1.4e} m')
print(f'Luminosity: L =  F*4*pi*r**2 = {L:1.4e} W')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
from scipy.constants import pi, parsec
f0 = 2.518021002E-8 # W/m/m
r = 10*parsec # m
L0 = 4*pi*r**2*f0 # W
print(f'r = 10 pc = {r:1.4e} m')
print(f'Conventional flux for zero apparent magnitude (m=0): F0={f0:1.4e} W/m/m')
print(f'Conventional luminosity for zero absolute magnitude (M=0): L0={L0:1.4e} W')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
from numpy import log10
from astropy.constants import L_sun, L_bol0

print(L_sun)
print(L_bol0)
L0 = 3.0128E28 # W
M = -2.5*log10(L_sun/L_bol0)
print(f'Absolute magnitude of the Sun: M_sun = {M:.2f}')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
from astropy.constants import L_bol0

M_sun = 4.74
print(L_bol0)
L = L_bol0*10**(-0.4*M_sun)
print(f'Luminosity of the Sun: L = {L:1.3e}')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: notes}
from numpy import log10
m=4.2
r=20 #pc
M = m + 5 - 5*log10(r)
print(f'M = m + 5 - 5*log10(r) = {m} + 5 - log10({r:.2f})={m}+5{-log10(r)} = {M}')
print(f'Absolute magnitude is {M:.2f}')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: skip}
#| label: fig-Rayleigh-Jeans-plot
#| fig-cap: "Spectral radiance plotted using the Rayleigh-Jeans law (@eq-rayleigh-jeans-law)."
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from astropy.constants import k_B, c, h
from matplotlib import pyplot as plt
from matplotlib import rc
rc('text',usetex=True)
rc('font',size=15)
rc('figure',figsize=(8,6), dpi=300)

l = np.linspace(0,3000,1000) # Wavelength, nm
T = 4000 # Temperature, Kelvin

B_rj=c.value* 2*k_B.value * T /(l*1e-9)**4; # 1e-9 to convert from nm to meters

plt.plot(l, B_rj, label="Rayleigh-Jeans law")
plt.xlabel(r'Wavelength $\lambda$ (nm)')
plt.ylabel(r'Spectral Radiance $B_\lambda(\lambda, T=4000K)$ (Wm$^{-3}$sr$^{-1}$)')
plt.title('Spectral radiance')
plt.ylim([0,1e15])
plt.grid()
plt.legend();
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: skip}
#| label: fig-Rayleigh-Wien-plot
#| fig-cap: "Comparison of the spectral radiance plotted using the Rayleigh-Jeans law (@eq-rayleigh-jeans-law) and the Wien's approximation (@eq-Wien_approx)."
B_wien=2 * h.value *c.value**2 * np.exp(-h.value*c.value/((l*1e-9)*k_B.value*T)) / (l*1e-9)**5;

plt.plot(l, B_rj, label="Rayleigh-Jeans law")
plt.plot(l, B_wien, label="Wien's approximation")

plt.xlabel(r'Wavelength $\lambda$ (nm)')
plt.ylabel(r'Spectral Radiance $B_\lambda(\lambda, T=4000K)$ (Wm$^{-3}$sr$^{-1}$)')
plt.title('Spectral radiance')
plt.ylim([0,5e12])
plt.grid()
plt.legend();
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: subslide}
#| label: fig-Rayleigh-Wien-Planck-plot
#| fig-cap: "Comparison of the spectral radiance plotted using the Rayleigh-Jeans law (@eq-rayleigh-jeans-law), Wien's approximation (@eq-Wien_approx) and Planck's law (@eq-planck-law)"
B_planck = 2*h.value*c.value**2 * (1.0/(np.exp(h.value*c.value/(k_B.value*T*(l*1e-9)))-1)) / (l*1e-9)**5;
plt.plot(l, B_planck, label="Planck's law")
plt.plot(l, B_wien, label="Wien's approximation")
plt.plot(l, B_rj, label="Rayleigh-Jeans law")

plt.xlabel(r'Wavelength $\lambda$ (nm)')
plt.ylabel(r'Spectral Radiance $B_\lambda(\lambda, T=4000K)$ (Wm$^{-3}$sr$^{-1}$)')
plt.title('Spectral radiance')
plt.ylim([0,5e12])
plt.grid()
plt.legend();
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: subslide}
#| label: fig-spectral-radiance-vs-temp
#| fig-cap: "Spectral radiance at different temperatures."

import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt

for T in [6000,5500,5000,4500,4000]:
    B_planck = 2*h.value*c.value**2 * (1.0/(np.exp(h.value*c.value/(k_B.value*T*(l*1e-9)))-1)) / (l*1e-9)**5
    plt.plot(l, B_planck, label='T='+str(T)+'K')


plt.xlabel(r'Wavelength $\lambda$ (nm)')
plt.ylabel(r'Spectral Radiance $B_\lambda(\lambda, T=4000K)$ (Wm$^{-3}$sr$^{-1}$)')
plt.title('Spectral radiance at different temperatures')
plt.grid()
plt.legend()
plt.show()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| slideshow: {slide_type: notes}
w = 0.0029 # constant in m K
lmax = 500e-9 # peak wavelength in metres
Tsun = w/lmax
print(f'The effective temperature of the Sun is approximately {Tsun:.3f} K')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
from numpy import log10
m = 3.54
r = 3.212 #pc
M = m + 5 - 5*log10(r)
print(f'Absolute magnitude: {M:.1f}')
#
#
#
#
#
#
#
#
#
#
#
#
L0 = 3.0128E28
L = L0*10**(-0.4*M)
print(f'Bolometric luminosity {L:.3e}')
#
#
#
#
#
#
#
#
#
#
#
#
c = 0.0029 # m K
wavel = 570E-9 # m 
T = c/wavel
print(f'Temperature: T={T:.3e}')
#
#
#
#
#
#
#
#
#
#
#
from numpy import sqrt, pi
sigma = 5.671E-8
R = sqrt(L/(4*pi*sigma*T**4))
print(f'Radius: R={R:.3e}')
#
#
#
#
#
#
#
#
#
#
#
