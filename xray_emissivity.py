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
    for X-ray emission calculations, assuming collisional ionization equilibrium. Uses scipy.interpolate.rbf for interpolation. 
    
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

        if self.rmf != None :
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

        # declare the Collisional Ionization Equilibrium cieion
        cie = pyatomdb.spectrum.CIESession()
        cie.set_broadening(True)

        # set the response (raw keyword tells pyatomdb it is not a real response file)
        if self.response == False: 
            print("Computing X-ray spectrum without instrumental response")
            cie.set_response(self.ebins, raw=True)

        else :

            try :
                path.exists(self.rmf)
            except OSError as error:
                logging.exception("RMF file "+str(self.rmf)+" not found!")
                raise error
            else :
                cie.set_response(self.rmf, self.arf)
                self.ebins = cie.ebins_out
                self.num_ebins = len(cie.ebins_out)
                self.emin = self.ebins[0]
                self.emax = self.ebins[-1]
        
        self.cie = cie

    def compute_cie_xray_spectrum (self, temperature, metallicity):

        Zlist = np.arange(31) #<- all the elements
        
        # This is the electron-electron bremstrahlung component alone
        #set all abundances to 1 (I need a full census of electrons in the plasma for e-e brems)
        self.cie.set_abund(Zlist[1:], metallicity)
        self.cie.set_eebrems(True)
        spec = self.cie.return_spectrum(temperature, dolines=False, docont=False, dopseudo=False)

        # turn off all the elements
        self.cie.set_abund(Zlist[1:], 0.0)
        # turn back on this element
        self.cie.set_abund(Zlist[1:], metallicity)
        # turn off e-e bremsstrahlung (avoid double counting)
        self.cie.set_eebrems(False)

        spec += self.cie.return_spectrum(temperature)
        
        #for Z in Zlist[1:]:
        #    self.cie.set_abund(Z, metallicity)
        #    spec += self.cie.return_spectrum(temperature)

        #print(spec)
        #print(spec)
        
        spectrum = spec
    
        return spectrum

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
        
        index = np.minimum(np.searchsorted(emarray, energy)-1, 13)
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
            Default is [0.001, 30.0]
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
        
        filename : str
            filename of the file. The saved file will be named "filename.hdf5"
            
        Returns
        -------
        
        None
        
        '''
        try:
            self.etable != None
            
        except Exception as exception:
            logging.exception("Emissivity table not set up yet!")
            raise exception

        else :
            if path.exists(filename+".hdf5") and overwrite == False:
                print("File already exists! Not saving file.")
            else :
                print("Saving emissivity table as "+filename+".hdf5")
                f = h5py.File(filename+".hdf5", "w")
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
            
        else :
            self.etable = f["emissivity_table"][:]
            self.tbins = f["temperature_array"][:]
            self.zbins = f["metallicity_array"][:]
            self.ebins = f["energy_array"][:]
            self.num_ebins = len(self.ebins)
            self.emin = self.ebins[0]
            self.emax = self.ebins[-1]

            f.close()

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
        
        if use_log10 :            
            lt = np.log10(temperature)
            lz = np.log10(metallicity)
        else :
            lt = (temperature)
            lz = (metallicity)

        if np.isscalar(lt) :   
            if lt < ltbins[0]:
                lt = ltbins[0]
            elif lt > ltbins[-1]:
                lt = ltbins[-1] 
            if lz < lzbins[0]:
                lz = lzbins[0]
            elif lz > lzbins[-1]:
                lz = lzbins[-1]
                
        else :
            
            shape = temperature.shape
            lt = lt.flatten()
            lz = lz.flatten()
            for i, val in enumerate(lt) :
                if val < ltbins[0]:
                    lt[i] = ltbins[0]
                if val > ltbins[-1]:
                    lt[i] = ltbins[-1]
            for i, val in enumerate(lz) :
                if val < lzbins[0]:
                    lz[i] = lzbins[0]
                if val > lzbins[-1]:
                    lz[i] = lzbins[-1]

                    
        xi = (lt, lz)
        result = interpn(points, values, xi, fill_value=0.0, method=method) 

        if np.isscalar(lt) :
            if result < 0:
                result = 0
        else: 
            #for i, val in enumerate(result) :
            #    if val < 0 :
            #        result[i] = 0.0
            
            result[result < 0] = 0.0
        
            result = np.reshape(result, shape)
                    
        return result