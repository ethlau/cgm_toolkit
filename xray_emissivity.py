import pyatomdb
import numpy
import h5py
import logging
from os import path
from scipy.interpolate import rbf

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
    def __init__(self, energy_range=[0.01,2.0], num_ebins=101, response=None, rmf=None, arf=None, use_physical_unit=True):
        '''
        Initialize the object
        
        Parameters
        ----------
        
        energy_range : list of two values
            sets the lower and upper bound of X-ray photon energy spectrum in keV (i.e. sets emin and emax)
        num_ebins : int
            sets the number of bins in X-ray photon energy spectrum
        response : string, optional
            whether to include instrumental response. Default is None
        rmf : string, optional if reponse is set to None
            Filename for the rmf reponse file, have to be set if response is set to True
        arf : string, optional
            Filename for the arf reponse file, optional
        use_physical_unit : boolean, optional
            True: use physical units for the emissivity: erg s^-1 cm^3
            False: use photon units for the emissivity: ph s^-1 cm^3
            Default is set to True
 
        '''
        self.emin = energy_range[0]
        self.emax = energy_range[1]
        self.num_ebins = num_ebins

        # create a set of energy bins (in keV) for the response. Note these are
        self.ebins = numpy.linspace(self.emin, self.emax, self.num_ebins)

        self.response = response
        self.rmf = rmf
        self.arf = arf

        self.etable = None
        self.tbins = None
        self.zbins = None
        self.num_tbins = None
        self.num_zbins = None
        self.temperature_range = None
        self.metallicity_range = None
        self.use_physical_unit = use_physical_unit

    def set_tbins(self, temperature_range, num_tbins, use_log10=True) :

        if use_log10 :

            temperature_min = numpy.log10(temperature_range[0])
            temperature_max = numpy.log10(temperature_range[1])
            tbins = numpy.logspace(temperature_min,temperature_max,num_tbins)
        else :

            temperature_min = (temperature_range[0])
            temperature_max = (temperature_range[1])
            tbins = numpy.linspace(temperature_min,temperature_max,num_tbins)

        self.num_tbins = num_tbins
        self.tbins = tbins

    def set_zbins(self, metallicity_range, num_zbins, use_log10=True) :

        if use_log10 :

            metallicity_min = numpy.log10(metallicity_range[0])
            metallicity_max = numpy.log10(metallicity_range[1])
            zbins = numpy.logspace(metallicity_min,metallicity_max,num_zbins)
        else :

            metallicity_min = (metallicity_range[0])
            metallicity_max = (metallicity_range[1])
            zbins = numpy.linspace(metallicity_min,metallicity_max,num_zbins)

        self.num_zbins = num_zbins
        self.zbins = zbins
    
        
    def compute_cie_xray_spectrum (self, temperature, metallicity):

        # declare the Collisional Ionization Equilibrium cieion
        cie = pyatomdb.spectrum.CIESession()
        cie.set_eebrems(True)
        cie.set_broadening(True)

        # set the response (raw keyword tells pyatomdb it is not a real response file)
        if self.response == None: 
            cie.set_response(self.ebins, raw=True)

        else :

            try :
                path.exists(rmf)
            except OSError as error:
                logging.exception("RMF file "+str(rmf)+" not found!")
                raise error
            else :
                cie.set_response(self.rmf, self.arf)

        kT = temperature # temperature in keV

        Zlist = numpy.arange(31) #<- all the elements
        cie.set_abund(Zlist[1:], metallicity)

        spectrum = cie.return_spectrum(kT)

        return spectrum

    def compute_specific_xray_emissivity (self, temperature, metallicity) :

        spectrum = self.compute_cie_xray_spectrum (temperature, metallicity)

        if self.use_physical_unit == True :
    
            #de = self.ebins[1:]-self.ebins[:-1]
            e =  0.5*(self.ebins[1:]+self.ebins[:-1])
            sum_spec = numpy.sum(e*spectrum)

        else :
            sum_spec = numpy.sum(spectrum)

        return sum_spec


    def tabulate_specific_xray_emissivity(self, temperature_range=[0.001, 20.0], metallicity_range=[0.02, 2.0], num_tbins=11, num_zbins=11) :

        '''
        Tabulate X-ray emissivity table
        
        Parameters
        ----------
        
        temperature_range : list, optional
            input temperature range in keV, list of two numbers: first marks the lower bound, last marks the upper bound. 
            Default is [0,001, 20.0]
        metallicity_range : list, optional
            input metalicity range in solar unit, list of two numbers: first marks the lower bound, last marks the upper bound.
            Default is [0.02, 2.0]
        num_tbins: int, optional
            number of bins temperature array, default is 11
        num_zbins: int, optional
            number of bins metallicity array, default is 11   
            
        Returns
        -------
        
        None
        
        '''
        
        self.temperature_range = temperature_range
        self.metallicity_range = metallicity_range
        self.set_tbins(temperature_range, num_tbins)
        self.set_zbins(metallicity_range, num_zbins)
        self.num_tbins = num_tbins
        self.num_zbins = num_zbins

        
        self.etable = numpy.zeros([self.num_tbins, self.num_zbins])

        for it, t in enumerate(self.tbins) :
            for iz, z in enumerate(self.zbins) :
                self.etable[it, iz] = self.compute_specific_xray_emissivity (t, z)

    def save_emissivity_table (self,filename) :
   
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
            if path.exists(filename+".hdf5") :
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

    def setup_interpolated_emissivity_table(self):

        '''
        Set up interpolated table from emissivity table
        '''
        
        try:
            self.etable != None
    
        except Exception as exception:
            logging.exception("Emissivity table not set up yet!")
            raise exception

        else :
        
            tbins2, zbins2 = numpy.meshgrid(numpy.log10(self.tbins), numpy.log10(self.zbins) , indexing='ij')
    
            self.e_rbf = rbf.Rbf(tbins2, zbins2, self.etable)


    def return_interpolated_emissivity (self, temperature, metallicity) :

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
            If response files are used, the unit is erg s^-1 cm^5 or ph s^-1 cm^5. 
        
        '''
        
        try:
            self.e_rbf != None

        except Exception as exception:
            logging.exception("Interpolated emissivity table not set up yet! Run `setup_interpolated_emissivity_table()` first!")
            raise exception

        else :
        
            return self.e_rbf(numpy.log10(temperature), numpy.log10(metallicity))
    

