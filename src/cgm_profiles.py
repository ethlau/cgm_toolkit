import numpy as np
import math

from astropy import units as un, constants
from astropy.cosmology import Planck18 as cosmo

from . import xray_emissivity
from . import ion_frac

Mpc_to_cm = 3.0857e24
kpc_to_cm = Mpc_to_cm/1.e3
pc_to_cm = Mpc_to_cm/1.e6
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
K_to_keV = 8.61733e-08 
Tcmb = 2.725 #K
v_c_ratio = 1.0e-3 # from Amodeo+21
h_planck = 6.6261e-27 # in cgs
kb = 1.3807e-16 # in cgs


class HaloProfile():

    def __init__ (self, redshift, radial_bin, pressure, density, metallicity, temperature=None) :

        self.redshift = redshift
        self.radial_bin = radial_bin # radius in kpc
        self.pressure = pressure # electron pressure in keV cm^-3
        self.density = density # electron number density in cm^-3
        self.metallicity = metallicity # metallicty in Zsolar
        #self.temperature = temperature # temperature in keV
        self.v_c_ratio = v_c_ratio

        if np.all(temperature) == True:
            self.temperature = temperature # temperature in keV
        else:
            self.temperature = pressure / density # pressure in keV cm^-3, density in cm^-3

    def set_v_c (self, v_c) :

        self.v_c_ratio = v_c
    

    def differential_volume (self, radius):
        
        '''
        Input : radius in kpc
        
        return volume of shell in kpc^3
        '''
        
        volume = np.zeros(radius.shape)
        
        volume[0] = (4.*math.pi/3.) * radius[0]**3

        for ir, r in enumerate(radius) :
            if ir > 0:
                volume[ir] = (4.*math.pi/3.) * (radius[ir]**3 - radius[ir-1]**3)
                
        return volume

    
    def abel_projection(self, r2D, prof) :
    
        '''
        return: np array with projected profile
        '''
        from scipy.interpolate import interp1d
        from scipy import integrate
    
        prof2D = np.zeros(len(r2D))
        prof3D = interp1d(r2D, prof, kind='linear', fill_value=None, bounds_error=False)

        for irp, rp in enumerate(r2D) :
            zmax = np.sqrt(max(r2D)**2 - rp**2)
            zbin_edge = np.linspace(0.0, zmax, 10000)
            zbin = 0.5*(zbin_edge[1:] + zbin_edge[:-1])
            dz = zbin_edge[1:]-zbin_edge[:-1]
            rprime = np.sqrt(rp**2+zbin**2)

            fz = 2*prof3D(rprime)
            integ = integrate.trapezoid(fz,zbin)
            prof2D[irp] = integ

        return prof2D 
    
    def spherical_dy_profile (self, radius) :

        '''
        Input: 
            radius: numpy array of radius in kpc. 
        
        return: spherical *differetial* compton-y profile in kpc^-2

        '''
        profile = np.interp(radius, self.radial_bin, self.pressure) * pe_factor

        profile *= sigma_T / m_e_keV * kpc_to_cm #unit = kpc^-1
        
        return profile

    def spherical_integrated_y_profile (self, radius) :

        '''
        Input: 
            radius: numpy array of radius in kpc. 
        
        return: spherical *integrated* compton-y profile in kpc^2

        '''
        
        dvol = self.differential_volume(radius) #kpc^3
        
        profile = self.spherical_dy_profile(radius) * dvol # kpc^-1 * kpc^3 -> kpc^2

        return profile

    def projected_y_profile (self, radius) :

        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: numpy array of y profile (dimensionless)

        '''
        
        sph_prof = self.spherical_dy_profile(radius) #kpc^-1

        profile = self.abel_projection(radius, sph_prof) 
        
        return profile

    def projected_tsz_profile (self, radius, frequency ) :

        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
            frequency: float in Hz. Frequency of microwave detector
        
        return: T_tsz in micro K

        '''


        x = h_planck * frequency / (kb * Tcmb)
       
        fnu = x * (1./np.tanh(x/2.0)) - 4.0

        profile = self.projected_y_profile(radius) #unitless

        profile *= fnu * Tcmb /  1.0e-6 
        
        return profile

        
    def spherical_xray_emissivity_profile (self, radius, etable='etable.hdf5') :

        '''
        
        input:
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: numpy array with spherical xray emissivity profile in erg/cm^3/s/sr

        '''
        xray = xray_emissivity.XrayEmissivity()
        xray.read_emissivity_table(etable)
        xcool = xray.return_interpolated_emissivity(self.temperature, self.metallicity)

        ne = self.density 
        nH = ne / 1.2
        em = xcool  * ne * nH / (1.+self.redshift)**4 / (4.0 * math.pi)

        profile = np.interp(radius, self.radial_bin, em)

        return profile
 
    def projected_xray_surface_brightness_profile (self, radius, etable='etable.hdf5') :

        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: numpy array with xray surface brightness profile in erg/cm^2/s/sr

        '''
        
        sph_prof = self.spherical_xray_emissivity_profile(radius, etable=etable) #erg/cm^3/s/sr

        profile = self.abel_projection(radius, sph_prof) * kpc_to_cm #erg/cm^2/s/sr

        return profile

    def projected_density_profile (self, radius) :

        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: numpy array of project number density profile (in cm^-2)

        '''
        
        sph_prof = self.density / pe_factor

        profile = self.abel_projection(radius, sph_prof) * kpc_to_cm
 
        return profile

    def projected_electron_density_profile (self, radius) :

        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: numpy array of project number density profile (in cm^-2)

        '''
        
        sph_prof = self.density
        
        profile = self.abel_projection(radius, sph_prof) * kpc_to_cm
 
        return profile

    def projected_ksz_profile (self, radius) :
        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: T_ksz in micro K

        '''


        profile = self.projected_electron_density_profile( radius )

        profile *= -Tcmb * sigma_T * v_c_ratio /1.0e-6

        return profile


    def return_ion_frac (self, temperature, ion, ion_table = './cgm_toolkit/data/ion_frac_CIE.txt' ):

        '''
        Input:
            temperature: temperture in keV
            ion: string name of ion, e.g., 'FeXXV'

        '''
        
        ifrac = ion_frac.IonFrac( ion_table = ion_table)
        
        temperature_in_K = temperature / K_to_keV

        frac = ifrac.return_ion_fraction (temperature_in_K, ion)
        
        return frac

    def return_ion_fraction_profile (self, radius, ion, ion_table = './cgm_toolkit/data/ion_frac.txt'):

        '''
        Input:
            temperature: temperture in keV
            ion: string name of ion, e.g., 'FeXXV'

        '''
        
        temperature = np.interp(radius, self.radial_bin, self.temperature)

        ion_frac = self.return_ion_frac (temperature, ion, ion_table=ion_table)

        return ion_frac



    def return_oxygen_column_densities (self, radius, ion, ion_table = './cgm_toolkit/data/ion_frac.txt'):

        '''
        Input:
            temperature: temperture in keV
            ion: string name of ion, e.g., 'FeXXV'

        '''

        a_o = 4.9e-4 # O/H solar
        
        temperature = np.interp(radius, self.radial_bin, self.temperature)

        ion_frac = self.return_ion_frac (temperature, ion, ion_table=ion_table)

        nH_profile = self.density / 1.2

        ion_density_profile = ion_frac * a_o * nH_profile *  self.metallicity 


        projected_profile = self.abel_projection(radius, ion_density_profile) * kpc_to_cm

        profile = projected_profile
        
        return profile

    def dispersion_measure_profile (self, radius) :

        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: numpy array of dispersion measure (in pc cm^-3)

        '''
        
        profile = self.projected_electron_density_profile(radius)

        profile *= 1./pc_to_cm

        return profile

        
