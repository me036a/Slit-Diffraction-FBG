import scipy
import numpy as np
from scipy.special import fresnel
import matplotlib
import matplotlib.pyplot as plt


def Fres_Diff(write_wavelength, slitwidth, x_pos, screen_pos):
    alpha_1 = (x_pos + slitwidth/2)*(np.sqrt(2/(write_wavelength*screen_pos)))
    alpha_2 = (x_pos - slitwidth/2)*(np.sqrt(2/(write_wavelength*screen_pos)))
    s1, c1 = fresnel(alpha_1) #Fresnel integrals
    s2, c2 = fresnel(alpha_2) #Fresnel integrals
    Intensity = 0.5*((c2-c1)**2+(s2-s1)**2)
    return Intensity

def Rouard1(wavelength, period, length, intensity):
    n_L = 1.5               #base refractive index 
    x_H = period/2           #length of half period
    x_L = x_H                 #xn and xp same length

    no_period = length/period
       
    ref_H = 0.0
    for i in range(0,np.int(no_period-1)):      
        
        n_H = n_L + 0.001*intensity[i]  #high refractive index is small change on base refractive index
        
        r_L = (n_H - n_L)/(n_H+n_L)
        r_H = (n_L - n_H)/(n_H+n_L)
        
        theta_L = (2*np.pi*n_L*x_L)/wavelength
        theta_H = (2*np.pi*n_H*x_H)/wavelength
        
        ref_L = (r_L + ref_H*np.exp(-2.0j * theta_L)) / (1+ r_L*ref_H*np.exp(-2.0j * theta_L))
        ref_H = (r_H + ref_L*np.exp(-2.0j * theta_H)) / (1+ r_H*ref_L*np.exp(-2.0j * theta_H))
        
    return np.abs(ref_H*np.conj(ref_H));  #returns reflection value from reflectivity

#

plt.close("all")

# Constants for Rouards Method

grating_period = 0.0005                                                         #grating period
wavelength = np.linspace(1.496e-3,1.504e-3,500)                                # wavelength range
d = 1.0                                                                         # slitwidth mm
width = np.int(3*d/grating_period)*grating_period                               # ensure that raneg in x direction is interger number of periods 
no_period = width/grating_period
x = np.linspace(-1.5*width/3,1.5*width/3,no_period)                                 # distance in x direction mm
length = np.max(x)- np.min(x)                                                   
lam = 2.66e-4       # wavelength (write) mm
#x = np.linspace(-1.5*d,1.5*d,5000) # distance in x direction mm
z = np.linspace(0,4,100) # log10 of the distance to the screen mm

I1 = np.zeros((len(x),len(z))) # set up array to receive the intensity data

# fix slitwidth, vary screen location

for i in range (0,len(z)-1):
    I1[:,i] = Fres_Diff(lam, d, x, 10**z[i])



zplot = 500 # distance at which you would  like to plot the intesity pattern mm

plt.figure(1)
plt.plot(x,I1[:,np.where(10**z>zplot-1)[0][0]])
plt.xlabel("x (mm)")
plt.ylabel("Intensity (au)")
plt.title("Intensity profile " +  str(int(10**z[np.where(10**z>zplot-1)[0][0]])) + "mm behind a slit of width " + str(d) + "mm")



Bragg_Spect = Rouard1(wavelength, grating_period, length, I1[:,np.where(10**z>zplot-1)[0][0]])

plt.figure(2)
plt.plot(wavelength,Bragg_Spect)

# plot limts
x_min = np.min(x)
x_max = np.max(x)
z_min = np.min(z)
z_max = np.max(z)


plt.figure(3)
plt.imshow(I1, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower', extent = [ 10**(z_min), 10**(z_max), x_min, x_max,] )
plt.xscale('log')
plt.xlabel('z (mm)')
plt.ylabel('x (mm)')

Bragg = np.zeros((len(wavelength),len(z)))

for i in range(0,len(z)-1):
    Bragg [:,i] = Rouard1(wavelength, grating_period, length, I1[:,i])
    print i



##
wave_min = np.min(wavelength)
wave_max = np.max(wavelength)

plt.figure(4)
plt.imshow(Bragg, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower', extent = [ 10**(z_min), 10**(z_max), wave_min, wave_max,] )
plt.xscale('log')
plt.xlabel('z (mm)')
plt.ylabel('Wavelength (microns)')


## fix screen position, vary slit width


d = np.linspace(0.2,1,100) # slitwidth  mm
x = np.linspace(-1.5*np.max(d),1.5*np.max(d),1000) # distance in x direction mm
z = 100 # distance to screen mm

I2 = np.zeros((len(x),len(d))) # set up array to receive the intensity data


for i in range (0,len(d)-1):
    I2[:,i] = Fres_Diff(lam, d[i], x, z)

# plot limits
x_min = np.min(x)
x_max = np.max(x)
d_min = np.min(d)
d_max = np.max(d)

plt.figure(5)
plt.imshow(I2, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower', extent = [ d_min, d_max, x_min, x_max,] )
plt.xlabel('Slitwidth (mm)')
plt.ylabel('x (mm)')
plt.title("Intensity as function of slitwidth" + "\n" + "assuming wavelength 266 nm and viewing distance " +  str(z) + "mm")

Bragg2 = np.zeros((len(wavelength),len(d)))

for i in range(0,len(d)-1):
    width = np.int(3*d[i]/grating_period)*grating_period                               # ensure that raneg in x direction is interger number of periods 
    no_period = width/grating_period
    x = np.linspace(-1.5*width/3,1.5*width/3,no_period)                                 # distance in x direction mm
    length = np.max(x)- np.min(x)
    I = Fres_Diff(lam, d[i], x, z)
    Bragg2[:,i] = Rouard1(wavelength, grating_period, length, I)
    print i


plt.figure(6)
plt.imshow(Bragg2, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower', extent = [ d_min, d_max, wave_min, wave_max,] )
plt.xlabel('Slitwidth (mm)')
plt.ylabel('Wavelength (microns)')
plt.title("Bragg spectrum as function of slitwidth" + "\n" + "assuming writing wavelength 266 nm and viewing distance " +  str(z) + "mm")

plt.show()