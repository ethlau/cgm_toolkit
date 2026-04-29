import numpy as np
import math

from astropy import units as un, constants
from astropy.cosmology import Planck18 as cosmo
from pathlib import Path

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
h_planck = 6.6261e-27 # in cgs
kb = 1.3807e-16 # in cgs


etable_file_path = Path(__file__).parent / "../data/etable_05_2keV.hdf5"
ion_frac_file_path = Path(__file__).parent / "../data/ion_frac_CIE.txt"


class HaloProfile():

    def __init__ (self, redshift, radial_bin, pressure, density, metallicity, temperature=None) :

        self.redshift = redshift
        self.radial_bin = radial_bin # radius in kpc
        self.pressure = pressure # electron pressure in keV cm^-3
        self.density = density # electron number density in cm^-3
        self.metallicity = metallicity # metallicty in Zsolar
        #self.temperature = temperature # temperature in keV

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

    
    def abel_projection(
        self,
        r,
        prof,
        rmax_factor=2.0,
        n_z=2048,
        fill_mode="powerlaw",
    ):
        """
        Project a spherical 3D profile along the line of sight.

        This version keeps the old interface:

            abel_projection(r, prof)

        but avoids forcing the projected profile to zero at the last bin by
        integrating out to rmax_factor * max(r).

        Parameters
        ----------
        r : array-like
            Radius array. Used both as projected radii and 3D profile radii.
        prof : array-like
            3D profile values.
        rmax_factor : float
        Outer integration radius in units of max(r).
        n_z : int
            Number of line-of-sight integration points.
        fill_mode : {"zero", "edge", "powerlaw"}
            How to handle radii beyond the tabulated range.

        Returns
        -------
        prof2D : ndarray
            Projected profile.
        """

        import numpy as np
        from scipy.interpolate import interp1d
        from scipy import integrate

        r = np.asarray(r, dtype=float)
        prof = np.asarray(prof, dtype=float)

        if np.any(r <= 0):
            raise ValueError("All radii must be positive.")

        if len(r) != len(prof):
            raise ValueError("r and prof must have the same length.")

        order = np.argsort(r)
        r = r[order]
        prof = prof[order]

        rmax_data = np.max(r)
        rmax = rmax_factor * rmax_data

        if rmax <= rmax_data:
            raise ValueError("rmax_factor must be greater than 1.")

        base_interp = interp1d(
            r,
            prof,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

        if fill_mode == "zero":

            def prof_interp(x):
                y = base_interp(x)
                return np.where(np.isfinite(y), y, 0.0)

        elif fill_mode == "edge":

            def prof_interp(x):
                y = base_interp(x)
                y = np.where(x < r[0], prof[0], y)
                y = np.where(x > r[-1], prof[-1], y)
                return y

        elif fill_mode == "powerlaw":
            # Estimate outer logarithmic slope from the last few radial bins.
            nfit = min(5, len(r))
            log_r_fit = np.log(r[-nfit:])
            log_p_fit = np.log(np.maximum(prof[-nfit:], 1e-300))

            slope, intercept = np.polyfit(log_r_fit, log_p_fit, 1)

            def prof_interp(x):
                x = np.asarray(x)
                y = base_interp(x)

                inner = x < r[0]
                outer = x > r[-1]

                y = np.where(inner, prof[0], y)

                outer_y = prof[-1] * (x / r[-1]) ** slope
                y = np.where(outer, outer_y, y)
    
                return y

        else:
            raise ValueError("fill_mode must be 'zero', 'edge', or 'powerlaw'.")

        prof2D = np.zeros_like(r)

        for i, R in enumerate(r):
            zmax = np.sqrt(rmax**2 - R**2)

            z = np.linspace(0.0, zmax, n_z)
            r_los = np.sqrt(R**2 + z**2)

            fz = 2.0 * prof_interp(r_los)
            prof2D[i] = integrate.trapezoid(fz, z)

        return prof2D
    
    
    def spherical_dy_profile (self, radius) :

        '''
        Input: 
            radius: numpy array of radius in kpc. 
        
        return: spherical *differetial* compton-y profile in kpc^-2

        '''
        profile = np.interp(radius, self.radial_bin, self.pressure)

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

        
    def spherical_xray_emissivity_profile (self, radius, etable=etable_file_path) :

        '''
        
        input:
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.

            etable: string of file name of tabulated emissivity table.
        
        return: numpy array with spherical xray emissivity profile
            
            unit depends on that used in the emissivity table. If in physical units, then it's erg/s/cm^3/steradian
            Otherwise it's in photons/s/cm^3/steradian

            If effective area of a given instrument was used in tabulating the emissivity table, then the units
            are either erg/s/cm/steradian, or photons/s/cm/steradian. 

        '''
        xray = xray_emissivity.XrayEmissivity()
        xray.read_emissivity_table(etable)
        xcool = xray.return_interpolated_emissivity(self.temperature, self.metallicity)

        ne = self.density 
        nH = ne / 1.2
        em = xcool  * ne * nH / (1.+self.redshift)**4 / (4.0 * math.pi)

        profile = np.interp(radius, self.radial_bin, em)

        return profile
 
    def projected_xray_surface_brightness_profile (self, radius, etable=etable_file_path) :

        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
        
        return: numpy array with xray surface brightness profile.
            unit depends on that used in the emissivity table. If in physical units, then it's erg/s/cm^2/steradian
            Otherwise it's in photons/s/cm^2/steradian

            If effective area of a given instrument was used in tabulating the emissivity table, then the units
            are either erg/s/steradian, or photons/s/steradian. 


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

    def projected_ksz_profile (self, radius, v_c_ratio = 1.e-3) :
        '''
        Input: 
            radius: numpy array of radius in kpc. Can be different from the one that initialize the class.
            v_c_ratio: ratio peculiar velocity to speed of light. Assume the same throughout the halo gas
        
        return: T_ksz in micro K

        '''


        profile = self.projected_electron_density_profile( radius )

        profile *= -Tcmb * sigma_T * v_c_ratio /1.0e-6

        return profile


    def return_ion_frac (self, temperature, ion, ion_table=ion_frac_file_path  ):

        '''
        Input:
            temperature: temperture in keV
            ion: string name of ion, e.g., 'FeXXV'

        '''
        
        ifrac = ion_frac.IonFrac( ion_table = ion_table)
        
        temperature_in_K = temperature / K_to_keV

        frac = ifrac.return_ion_fraction (temperature_in_K, ion)
        
        return frac

    def return_ion_fraction_profile (self, radius, ion, ion_table=ion_frac_file_path ):

        '''
        Input:
            temperature: temperture in keV
            ion: string name of ion, e.g., 'FeXXV'

        '''
        
        temperature = np.interp(radius, self.radial_bin, self.temperature)

        ion_frac = self.return_ion_frac (temperature, ion, ion_table=ion_table)

        return ion_frac



    def return_oxygen_column_densities (self, radius, ion, ion_table=ion_frac_file_path ):

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

        
