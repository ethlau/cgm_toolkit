import numpy as np
import math

import xray_emissivity

Mpc_to_cm = 3.0857e24
sigma_T = 6.6524587158e-25 # in cm^2
m_e_keV = 510.9989461 # in keV
clight = 29979245800.0 #in cm/s
m_e = 9.11e-28 # electron mass in g
m_p = 1.6726e-24 # in g
Msun = 1.98892e33
mu = 0.58824; # X=0.76 assumed
mu_e = 1.136; # X=0.76 assumed
XH = 0.76 #primordial hydrogen fraction
pe_factor = (2.0*XH+2.0)/(5.0*XH+3.0) #/conversion factor from gas pressure to electron pressure

class HaloProfile():

    def __init__ (self, mass, redshift, radial_bin, pressure, density, metallicity, temperature=None) :

        self.mass = mass
        self.redshift = redshift
        self.radial_bin = radial_bin
        self.pressure = pressure
        self.density = density
        self.metallicity = metallicity

        if temperature:
            self.temperature = temperature # temperature in keV
        else :
            self.temperature = pressure / density # pressure in keV cm^-3, density in cm^-3

    
    def spherical_y_profile (self, radius) :

        '''
        return: spherical compton-y profile in Mpc^-1

        '''
        profile = np.interp(radius, self.radial_bin, self.pressure)

        profile *= pe_factor * sigma_T / m_e_keV * Mpc_to_cm  #unit = Mpc^-1

        return profile

        
    def spherical_xray_emissivity_profile (self, radius, etable='etable_erosita.hdf5') :

        '''
        return: spherical xray emissivity profile in erg/cm^3/s/sr

        '''
        xray = xray_emissivity.XrayEmissivity()
        xray.read_emissivity_table(etable)
        xcool = xray.return_interpolated_emissivity(self.temperature, self.metallicity)

        ne = self.density * mu / mu_e
        nH = ne / 1.2
        em = xcool  * ne * nH / (1.+self.redshift)**4 / (4.0 * math.pi)

        profile = np.interp(radius, self.radial_bin, em)

        return profile
