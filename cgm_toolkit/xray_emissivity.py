import pyatomdb
import numpy as np
import h5py
import logging
from os import path
from scipy.interpolate import interpn

from tqdm import tqdm

class XrayEmissivity():

    '''
    
    Class for tabulating X-ray emissivity values from temperature (in keV) and metallcity (in solar units). Uses pyatomdb module 
    for X-ray emission calculations, assuming collisional ionization equilibrium. Uses scipy.interpolate.interpn for interpolation. 
    
    Attributes
    ----------
    
    emin : double
        lower bound of X-ray photon energy spectrum in keV
    emax : double
        upper bound of X-ray photon energy spectrum in keV
    num_ebins : int
        number of bins in the X-ray photon energy spectrum
    ebins : numpy array
        array of X-ray photon energy spectrum in keV 
        
    response : boolean
        whether to include instrumental response. Default is False
    rmf : string
        Filename for the rmf reponse file, have to be set if response is set to True
    arf : string
        Filename for the arf reponse file, optional. If set, the units of the 
        
    etable : numpy array
        2D numpy array to store tabulated emissivity values. Its shape is determined by tbins and zbins.
    tbins : numpy array
        1D numpy array of the input values of temperatures in keV for emissivity tabulation
    zbins : numpy array
        1D numpy array of the input values of metallicity (in solar unit) for emissivity tabulation

    Methods
    -------
    tabulate_specific_xray_emissivity()
    save_emissivity_table (filename)
    read_emissivity_table (filename)
    setup_interpolated_emissivity_table()
    return_interpolated_emissivity (temperature, metallicity)
    
    
    '''
    def __init__(self, energy_range=[0.01,2.0], num_ebins=101, rmf=None, arf=None, use_energy_unit=True):
        '''
        Initialize the object
        
        Parameters
        ---------- 
        energy_range : list of two values
            sets the lower and upper bound of X-ray photon energy spectrum in keV (i.e. sets emin and emax)
        num_ebins : int
            sets the number of bins in X-ray photon energy spectrum
        rmf : string, optional if reponse is set to None
            Filename for the rmf reponse file, have to be set if response is set to True
        arf : string, optional
            Filename for the arf reponse file, optional
        use_energy_unit : boolean, optional
            True: use energy units for the emissivity: erg s^-1 cm^3
            False: use photon units for the emissivity: ph s^-1 cm^3
            Default is set to True
 
        '''
        self.emin = energy_range[0]
        self.emax = energy_range[1]
        self.num_ebins = num_ebins

        # create a set of energy bins (in keV) for the response. Note these are
        self.ebins = np.linspace(self.emin, self.emax, self.num_ebins)

        self.response = False
        self.rmf = rmf
        self.arf = arf

        # FIX: use 'is not None' instead of '!= None' to avoid ambiguous
        # truth-value comparison (would raise ValueError on numpy arrays).
        if self.rmf is not None:
            self.response = True
        
        self.etable = None
        self.tbins = None
        self.zbins = None
        self.num_tbins = None
        self.num_zbins = None
        self.temperature_range = None
        self.metallicity_range = None
        self.use_energy_unit = use_energy_unit
        
        self.cie = None

    def set_tbins(self, temperature_range, num_tbins, use_log10=True) :

        if use_log10 :

            temperature_min = np.log10(temperature_range[0])
            temperature_max = np.log10(temperature_range[1])
            tbins = np.logspace(temperature_min,temperature_max,num_tbins)
        else :

            temperature_min = (temperature_range[0])
            temperature_max = (temperature_range[1])
            tbins = np.linspace(temperature_min,temperature_max,num_tbins)

        self.num_tbins = num_tbins
        self.tbins = tbins

    def set_zbins(self, metallicity_range, num_zbins, use_log10=True) :

        if use_log10 :

            metallicity_min = np.log10(metallicity_range[0])
            metallicity_max = np.log10(metallicity_range[1])
            zbins = np.logspace(metallicity_min,metallicity_max,num_zbins)
        else :

            metallicity_min = (metallicity_range[0])
            metallicity_max = (metallicity_range[1])
            zbins = np.linspace(metallicity_min,metallicity_max,num_zbins)

        self.num_zbins = num_zbins
        self.zbins = zbins
    
        
    def setup_cie_xray_spectrum (self):

        # declare the Collisional Ionization Equilibrium session
        cie = pyatomdb.spectrum.CIESession()
        cie.set_broadening(True)

        # set the response (raw keyword tells pyatomdb it is not a real response file)
        if self.response == False: 
            print("Computing X-ray spectrum without instrumental response")
            cie.set_response(self.ebins, raw=True)

        else :

            # path.exists returns False rather than raising, so check it explicitly
            if not path.exists(self.rmf):
                raise FileNotFoundError("RMF file "+str(self.rmf)+" not found!")
            if self.arf is not None and not path.exists(self.arf):
                raise FileNotFoundError("ARF file "+str(self.arf)+" not found!")

            # pyatomdb expects arf=False when there is no ARF
            arf = self.arf if self.arf is not None else False
            cie.set_response(self.rmf, arf)
            self.ebins = cie.ebins_out
            self.num_ebins = len(cie.ebins_out)
            self.emin = self.ebins[0]
            self.emax = self.ebins[-1]
        
        self.cie = cie

    def compute_cie_xray_spectrum (self, temperature, metallicity):

        Zlist = np.arange(31)  # all the elements

        # Metallicity scales only the metals (Z >= 3). H and He stay at their
        # solar values, otherwise the H/He continuum (which dominates at low
        # metallicity and high temperature) would wrongly scale with Z.
        self.cie.set_abund(Zlist[1:3], 1.0)
        self.cie.set_abund(Zlist[3:], metallicity)

        # One call gives lines + continuum + pseudo-continuum + e-e
        # bremsstrahlung, with the e-e brems electron census taken from the
        # abundances above. No double counting happens within a single call.
        self.cie.dolines    = True
        self.cie.docont     = True
        self.cie.dopseudo   = True
        self.cie.do_eebrems = True
        spec = self.cie.return_spectrum(temperature)

        return spec

    def wabs (self, energy) :
        #wabs galactic absorption model */
        # E in keV

        emarray = np.array([0.1, 0.284, 0.4, 0.532, 0.707, 0.867, 1.303, 1.840,
                          2.471, 3.210, 4.038, 7.111, 8.331, 10.0])

        c0 = np.array([17.3, 34.6, 78.1, 71.4, 95.5, 308.9, 120.6, 141.3,
                         202.7, 342.7, 352.2, 433.9, 629.0, 701.2])

        c1 = np.array([608.1, 267.9, 18.8, 66.8, 145.8, -380.6, 169.3,
                        146.8, 104.7, 18.7, 18.7, -2.4, 30.9, 25.2])

        c2 = np.array([-2150., -476.1 ,4.3, -51.4, -61.1, 294.0, -47.7,
                        -31.5, -17.0, 0.0, 0.0, 0.75, 0.0, 0.0])
        
        # clip to [0, 13]: energies below 0.1 keV would otherwise give index -1,
        # which wraps around to the coefficients of the 10 keV segment
        index = np.clip(np.searchsorted(emarray, energy)-1, 0, 13)
        sigma = (c0[index]+c1[index]*energy+c2[index]*energy*energy)*1.0e-24/energy**3
        return sigma
    
    def compute_xray_emissivity (self, temperature, metallicity, nH) :

        spectrum = self.compute_cie_xray_spectrum (temperature, metallicity)

        if nH > 0:
            e =  0.5*(self.ebins[1:]+self.ebins[:-1])
            absorption = np.exp( - self.wabs(e) * nH*1e20 );
            spectrum *= absorption
        
        if self.use_energy_unit == True :
            e =  0.5*(self.ebins[1:]+self.ebins[:-1])
            sum_spec = np.sum(e*spectrum*pyatomdb.const.ERG_KEV)
        else :
            sum_spec = np.sum(spectrum)

        return sum_spec


    def tabulate_xray_emissivity(self, temperature_range=[0.001, 100.0], metallicity_range=[0.01, 10.0], nH = 0.0, num_tbins=51, num_zbins=11, use_log10=True) :

        '''
        Tabulate X-ray emissivity table
        
        Parameters
        ----------
        temperature_range : list, optional
            input temperature range in keV, list of two numbers: first marks the lower bound, last marks the upper bound. 
            Default is [0.001, 100.0]
        metallicity_range : list, optional
            input metalicity range in solar unit, list of two numbers: first marks the lower bound, last marks the upper bound.
            Default is [0.01, 10.0]
        num_tbins: int, optional
            number of bins temperature array, default is 51
        num_zbins: int, optional
            number of bins metallicity array, default is 11   
        nH : double, optional
            Hydrogen column density in 1e20 cm^-2 for Galactic absorption 

        Returns
        ------- 
        None
        
        '''
        
        self.temperature_range = temperature_range
        self.metallicity_range = metallicity_range
        self.set_tbins(temperature_range, num_tbins, use_log10=use_log10)
        self.set_zbins(metallicity_range, num_zbins, use_log10=use_log10)
        self.num_tbins = num_tbins
        self.num_zbins = num_zbins

        self.etable = np.zeros([self.num_tbins, self.num_zbins])

        self.setup_cie_xray_spectrum()
        
        for it, t in enumerate(tqdm(self.tbins)) :
            for iz, z in enumerate(self.zbins):
                self.etable[it, iz] = self.compute_xray_emissivity (t, z, nH)

    def save_emissivity_table (self,filename, overwrite=True) :
   
        '''
        Save X-ray emissivity table to file in hdf5 format
        
        Parameters
        ----------
        filename : string
            filename of the file. The saved file will be named "filename.hdf5"
            
        Returns
        -------        
        None
        
        '''
        # FIX: 'self.etable != None' raises ValueError when etable is a numpy
        # array because numpy performs element-wise comparison and the result
        # is an array, which has no single truth value.  Use 'is not None'.
        if self.etable is None:
            raise RuntimeError("Emissivity table not set up yet!")

        if (not filename.endswith('.hdf5')):
            filename = filename + '.hdf5'
        if path.exists(filename) and overwrite == False:
            print("File already exists! Not saving file.")
        else :
            print("Saving emissivity table as "+filename)
            f = h5py.File(filename, "w")
            dset1 = f.create_dataset("emissivity_table", self.etable.shape, dtype='f', data=self.etable)
            dset2 = f.create_dataset("temperature_array", self.tbins.shape, dtype='f', data=self.tbins)
            dset3 = f.create_dataset("metallicity_array", self.zbins.shape, dtype='f', data=self.zbins)
            dset4 = f.create_dataset("energy_array", self.ebins.shape, dtype='f', data=self.ebins)      
            f.close()

    def read_emissivity_table (self,filename) :

        
        '''
        Read previously saved X-ray emissivity table to file in hdf5 format
        
        Parameters
        ----------
        
        filename : str
            filename of the file. The saved file should be named as "filename.hdf5"
            
        Returns
        -------
        
        None
        
        '''
        
        try:
            f = h5py.File(filename, 'r')
        except OSError as error:
            logging.exception("File "+str(filename)+" not found!")
            raise error

        with f:
            self.etable = f["emissivity_table"][:]
            self.tbins = f["temperature_array"][:]
            self.zbins = f["metallicity_array"][:]
            self.ebins = f["energy_array"][:]
            self.num_ebins = len(self.ebins)
            self.emin = self.ebins[0]
            self.emax = self.ebins[-1]

    def return_interpolated_emissivity (self, temperature, metallicity, method='splinef2d', use_log10=True) :

        '''
        Returns emissivity value from interpolated emissivity table
        
        Parameters
        ----------
        temperature : double
            temperature in keV
        metallicity: double
            metallicity in solar units
        
        Returns
        -------
        emissivity : double
            emissivity in either erg s^-1 cm^3 if use_physical_unit has been set to True (default), or ph s^-1 cm^3 otherwise.
            If response files are used, the unit is erg s^-1 cm^5 in physical unit or ph s^-1 cm^5 otherwise. 
        
        '''
        
        if use_log10 :
            ltbins = np.log10(self.tbins)
            lzbins = np.log10(self.zbins)
        else :
            ltbins = (self.tbins)
            lzbins = (self.zbins)    
        points = (ltbins, lzbins)
        values = self.etable
        
        scalar_input = np.isscalar(temperature) and np.isscalar(metallicity)

        # accept lists, and allow e.g. an array of temperatures with one metallicity
        temperature, metallicity = np.broadcast_arrays(
            np.asarray(temperature, dtype=float), np.asarray(metallicity, dtype=float))
        shape = temperature.shape

        if use_log10 :
            lt = np.log10(temperature)
            lz = np.log10(metallicity)
        else :
            lt = temperature
            lz = metallicity

        # clamp to the table range (NaN passes through and gives fill_value)
        lt = np.clip(lt.ravel(), ltbins[0], ltbins[-1])
        lz = np.clip(lz.ravel(), lzbins[0], lzbins[-1])

        xi = (lt, lz)
        result = interpn(points, values, xi, fill_value=0.0, method=method, bounds_error=False)
        result[result < 0] = 0.0

        if scalar_input :
            return result[0]
        return np.reshape(result, shape)
