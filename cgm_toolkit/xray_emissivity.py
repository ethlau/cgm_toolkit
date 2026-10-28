import errno
import logging
import os
from pathlib import Path

import h5py
import numpy as np
import pyatomdb
from scipy.interpolate import interpn
from tqdm import tqdm

tbabs_data_path = Path(__file__).parent / "../data/xsect_tbabs_wilm.npz"


angr_abund = np.array(
    [
        1.0,
        0.0977,
        1.45e-11,
        1.41e-11,
        3.98e-10,
        3.63e-4,
        1.12e-4,
        8.51e-4,
        3.63e-8,
        1.23e-4,
        2.14e-6,
        3.80e-5,
        2.95e-6,
        3.55e-5,
        2.82e-7,
        1.62e-5,
        3.16e-7,
        3.63e-6,
        1.32e-7,
        2.29e-6,
        1.26e-9,
        9.97e-8,
        1.00e-8,
        4.68e-7,
        2.45e-7,
        4.68e-5,
        8.32e-8,
        1.78e-6,
        1.62e-8,
        3.98e-8,
    ]
)

apsl_abund = np.array(
    [
        1.0,
        0.0851,
        1.12e-11,
        2.40e-11,
        5.01e-10,
        2.69e-4,
        6.76e-5,
        4.90e-4,
        3.63e-8,
        8.51e-5,
        1.74e-6,
        3.98e-5,
        2.82e-6,
        3.24e-5,
        2.57e-7,
        1.32e-5,
        3.16e-7,
        2.51e-6,
        1.07e-7,
        2.19e-6,
        1.41e-9,
        8.91e-8,
        8.51e-9,
        4.37e-7,
        2.69e-7,
        3.16e-5,
        9.77e-8,
        1.66e-6,
        1.55e-8,
        3.63e-8,
    ]
)


class XrayEmissivity:
    """
    Class for tabulating X-ray emissivity values as a function of temperature
    (keV), redshift, and metallicity (solar units). Uses pyatomdb for CIE
    X-ray emission calculations and scipy.interpolate.interpn for interpolation.

    The emissivity table is 3-dimensional with axes:
        [num_tbins, num_rbins, num_mbins]

    Attributes
    ----------
    emin : float
        Lower bound of X-ray photon energy spectrum in keV (observer frame).
    emax : float
        Upper bound of X-ray photon energy spectrum in keV (observer frame).
    num_ebins : int
        Number of bins in the X-ray photon energy spectrum.
    ebins : ndarray
        Energy bin edges in keV (observer frame).
    response : bool
        Whether to include an instrumental response. Default False.
    rmf : str or None
        Filename for the RMF response file. Required when response is True.
    arf : str or None
        Filename for the ARF response file. Optional.
    etable : ndarray, shape (num_tbins, num_rbins, num_mbins)
        Tabulated emissivity values.
    tbins : ndarray
        Temperature bin centres in keV.
    rbins : ndarray
        Redshift bin centres (1 + z).
    mbins : ndarray
        Metallicity bin centres in solar units.

    Methods
    -------
    tabulate_xray_emissivity(...)
    save_emissivity_table(filename)
    read_emissivity_table(filename)
    return_interpolated_emissivity(temperature, redshift, metallicity)
    """

    def __init__(self):
        self.emin = None
        self.emax = None
        self.num_ebins = None
        self.ebins = None

        self.rmf = None
        self.arf = None
        self.response = False

        self.etable = None
        self.tbins = None
        self.rbins = None  # stored as (1 + z) for log-spacing convenience
        self.mbins = None
        self.num_tbins = None
        self.num_rbins = None
        self.num_mbins = None
        self.temperature_range = None
        self.redshift_range = None
        self.metallicity_range = None
        self.use_energy_unit = True

        self.cie = None

    # ---------------------------------------------------------------------- #
    # Bin setup
    # ---------------------------------------------------------------------- #

    def set_tbins(self, temperature_range, num_tbins, use_log10=True):
        """Set temperature bin array in keV."""
        if use_log10:
            tbins = np.logspace(
                np.log10(temperature_range[0]),
                np.log10(temperature_range[1]),
                num_tbins,
            )
        else:
            tbins = np.linspace(temperature_range[0], temperature_range[1], num_tbins)
        self.num_tbins = num_tbins
        self.tbins = tbins

    def set_rbins(self, redshift_range, num_rbins, use_log10=True):
        """
        Set redshift bin array. Bins are stored as (1 + z) so that
        log-spacing remains well-defined at z = 0.
        """
        zmin = 1.0 + redshift_range[0]
        zmax = 1.0 + redshift_range[1]
        if use_log10:
            rbins = np.logspace(np.log10(zmin), np.log10(zmax), num_rbins)
        else:
            rbins = np.linspace(zmin, zmax, num_rbins)
        self.num_rbins = num_rbins
        self.rbins = rbins

    def set_mbins(self, metallicity_range, num_mbins, use_log10=True):
        """Set metallicity bin array in solar units."""
        if use_log10:
            mbins = np.logspace(
                np.log10(metallicity_range[0]),
                np.log10(metallicity_range[1]),
                num_mbins,
            )
        else:
            mbins = np.linspace(metallicity_range[0], metallicity_range[1], num_mbins)
        self.num_mbins = num_mbins
        self.mbins = mbins

    # ---------------------------------------------------------------------- #
    # CIE session
    # ---------------------------------------------------------------------- #

    def setup_cie_xray_spectrum(self, redshift, nH=0.0):
        """
        Initialise a CIESession for a given redshift.

        The spectrum is always computed on observer-frame bins: pyatomdb
        evaluates the model at the rest-frame energies E_obs * (1 + z), so
        each bin holds the rest-frame photons that land in that observed bin.

        Parameters
        ----------
        redshift : float
            Source redshift z.
        nH : float, optional
            Galactic hydrogen column density in cm^-2. With a response, the
            absorption is folded into the ARF so it acts on the incoming
            (observer-frame) photon energy, before redistribution by the RMF.
        """
        cie = pyatomdb.spectrum.CIESession()
        cie.set_broadening(True)

        if not self.response:
            print("Computing X-ray spectrum without instrumental response")
            cie.set_response(self.ebins * (1.0 + redshift), raw=True)
        else:
            # Path.exists() returns False rather than raising, so check explicitly
            for fname in (self.rmf, self.arf):
                if fname is not None and not Path(fname).exists():
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), fname
                    )
            # pyatomdb expects arf=False when there is no ARF
            cie.set_response(self.rmf, self.arf if self.arf is not None else False)

            # Model input bins from the RMF are observer-frame energies.
            # Evaluate the model at the matching rest-frame energies so the
            # ARF/RMF act on the observed photons.
            specbins_obs = np.array(cie.specbins, dtype=float)
            cie.specbins = specbins_obs * (1.0 + redshift)

            if nH > 0:
                e_in = 0.5 * (specbins_obs[1:] + specbins_obs[:-1])
                cie.arf = cie.arf * np.exp(-self.tbabs(e_in) * nH)

            # Output channels (observer frame) and the channels inside the band
            self.channel_ebins = np.array(cie.ebins_out, dtype=float)
            chan_mid = 0.5 * (self.channel_ebins[1:] + self.channel_ebins[:-1])
            self.channel_mask = (chan_mid >= self.emin) & (chan_mid <= self.emax)

        self.cie = cie

    # ---------------------------------------------------------------------- #
    # Spectrum and emissivity computation
    # ---------------------------------------------------------------------- #

    def compute_cie_xray_spectrum(self, temperature, metallicity, redshift):
        """
        Compute the CIE X-ray spectrum for given plasma conditions.

        H and He abundances are scaled by the APSL/ANGR ratio only;
        heavier elements are additionally scaled by metallicity.

        Parameters
        ----------
        temperature : float
            Plasma temperature in keV.
        metallicity : float
            Metallicity in solar units.
        redshift : float
            Source redshift z. Unused here: the redshift is applied to the
            energy grid in :meth:`setup_cie_xray_spectrum`. Kept for API
            compatibility.

        Returns
        -------
        spec : ndarray
            Rest-frame photon emissivity per observer-frame bin, in
            ph cm^3 s^-1 bin^-1 (raw), or counts cm^5 s^-1 channel^-1 with a
            response.
        """
        Zlist = np.arange(31)

        # Abundance ratios: APSL values relative to the Anders & Grevesse
        # reference built into pyatomdb.
        # H and He are not scaled by metallicity; all heavier elements are.
        self.cie.set_abund(Zlist[1:3], apsl_abund[0:2] / angr_abund[0:2])
        self.cie.set_abund(
            Zlist[3:], metallicity * apsl_abund[2:] / angr_abund[2:]
        )

        # One call gives lines + continuum + pseudo-continuum + e-e
        # bremsstrahlung (electron census from the abundances above); there
        # is no double counting within a single call.
        self.cie.dolines    = True
        self.cie.docont     = True
        self.cie.dopseudo   = True
        self.cie.do_eebrems = True
        spec = self.cie.return_spectrum(temperature)

        return spec

    def wabs(self, energy):
        """
        Wisconsin (Morrison & McCammon 1983) photoelectric absorption.

        Parameters
        ----------
        energy : array-like
            Photon energies in keV.

        Returns
        -------
        sigma : ndarray
            Cross-section in cm^2 per H atom.
        """
        emarray = np.array(
            [0.1, 0.284, 0.4, 0.532, 0.707, 0.867, 1.303, 1.840,
             2.471, 3.210, 4.038, 7.111, 8.331, 10.0]
        )
        c0 = np.array(
            [17.3, 34.6, 78.1, 71.4, 95.5, 308.9, 120.6, 141.3,
             202.7, 342.7, 352.2, 433.9, 629.0, 701.2]
        )
        c1 = np.array(
            [608.1, 267.9, 18.8, 66.8, 145.8, -380.6, 169.3, 146.8,
             104.7, 18.7, 18.7, -2.4, 30.9, 25.2]
        )
        c2 = np.array(
            [-2150.0, -476.1, 4.3, -51.4, -61.1, 294.0, -47.7, -31.5,
             -17.0, 0.0, 0.0, 0.75, 0.0, 0.0]
        )
        # clip to [0, 13]: below 0.1 keV the index would be -1, which wraps
        # around to the 10 keV coefficients
        index = np.clip(np.searchsorted(emarray, energy) - 1, 0, 13)
        return (
            (c0[index] + c1[index] * energy + c2[index] * energy**2)
            * 1.0e-24
            / energy**3
        )

    def tbabs(self, energy, tbabs_data_path=tbabs_data_path):
        """
        Tuebingen-Boulder ISM absorption model (Wilms et al. 2000).
        Data extracted from jaxspec (https://github.com/renecotyfanboy/jaxspec).

        Parameters
        ----------
        energy : array-like
            Photon energies in keV.

        Returns
        -------
        sigma : ndarray
            Cross-section in cm^2 per H atom.
        """
        tbabs_data = np.load(tbabs_data_path)
        sigma = np.interp(
            energy, tbabs_data["energy"], tbabs_data["sigma"], left=1e9, right=0.0
        )
        return sigma * 1e-22

    def compute_xray_emissivity(self, temperature, metallicity, redshift, nH):
        """
        Compute the band-integrated X-ray emissivity for a single set of
        plasma conditions.

        Parameters
        ----------
        temperature : float
            Plasma temperature in keV.
        metallicity : float
            Metallicity in solar units.
        redshift : float
            Source redshift z.
        nH : float
            Hydrogen column density in cm^-2 for Galactic absorption.

        Returns
        -------
        emissivity : float
            Emissivity in erg s^-1 cm^3 (use_energy_unit=True) or
            ph s^-1 cm^3 (use_energy_unit=False), with photon energies
            measured in the observer frame.
        """
        spectrum = self.compute_cie_xray_spectrum(temperature, metallicity, redshift)

        if self.response:
            # Channels inside the observed band. Absorption was already
            # folded into the ARF in setup_cie_xray_spectrum.
            e = 0.5 * (self.channel_ebins[1:] + self.channel_ebins[:-1])
            e = e[self.channel_mask]
            spectrum = spectrum[self.channel_mask]
        else:
            e = 0.5 * (self.ebins[1:] + self.ebins[:-1])
            if nH > 0:
                spectrum = spectrum * np.exp(-self.tbabs(e) * nH)

        # Photons are weighted by their observer-frame energy. The table is
        # therefore (rest-frame photon rate) x E_obs, and surface brightness
        # must be dimmed by (1+z)^3, not (1+z)^4 (see diffgas.xray.xray_dimming).
        if self.use_energy_unit:
            return np.sum(e * spectrum * pyatomdb.const.ERG_KEV)
        else:
            return np.sum(spectrum)

    # ---------------------------------------------------------------------- #
    # Tabulation
    # ---------------------------------------------------------------------- #

    def tabulate_xray_emissivity(
        self,
        energy_range=None,
        num_ebins=101,
        rmf=None,
        arf=None,
        use_energy_unit=True,
        temperature_range=None,
        redshift_range=None,
        metallicity_range=None,
        nH=0.0,
        num_tbins=51,
        num_rbins=11,
        num_mbins=11,
        use_log10=True,
    ):
        """
        Build the 3D emissivity table over temperature, redshift, and
        metallicity.

        Parameters
        ----------
        energy_range : list of two floats, optional
            Lower and upper bounds of the observer-frame energy band in keV.
            Default [0.5, 2.0].
        num_ebins : int, optional
            Number of energy bins. Default 101.
        rmf : str or None, optional
            RMF response filename. If provided, instrumental response is applied.
        arf : str or None, optional
            ARF response filename. Optional companion to rmf.
        use_energy_unit : bool, optional
            True → erg s^-1 cm^3; False → ph s^-1 cm^3. Default True.
        temperature_range : list of two floats, optional
            Min and max plasma temperature in keV. Default [0.0009, 30.0].
        redshift_range : list of two floats, optional
            Min and max source redshift. Default [0.0, 3.0].
        metallicity_range : list of two floats, optional
            Min and max metallicity in solar units. Default [0.01, 10.0].
        nH : float, optional
            Hydrogen column density in cm^-2 for Galactic tbabs absorption.
            Default 0 (no absorption).
        num_tbins : int, optional
            Number of temperature bins. Default 51.
        num_rbins : int, optional
            Number of redshift bins. Default 11.
        num_mbins : int, optional
            Number of metallicity bins. Default 11.
        use_log10 : bool, optional
            Logarithmic spacing for all three axes. Default True.

        Returns
        -------
        None
        """
        if energy_range is None:
            energy_range = [0.5, 2.0]
        if temperature_range is None:
            temperature_range = [0.0009, 30.0]
        if redshift_range is None:
            redshift_range = [0.0, 3.0]
        if metallicity_range is None:
            metallicity_range = [0.01, 10.0]

        # Energy grid (observer frame)
        self.emin = energy_range[0]
        self.emax = energy_range[1]
        self.num_ebins = num_ebins
        self.ebins = np.linspace(self.emin, self.emax, self.num_ebins)

        # Response
        self.rmf = rmf
        self.arf = arf
        self.response = self.rmf is not None

        self.use_energy_unit = use_energy_unit

        # Parameter grids
        self.temperature_range = temperature_range
        self.redshift_range = redshift_range
        self.metallicity_range = metallicity_range
        self.set_tbins(temperature_range, num_tbins, use_log10=use_log10)
        self.set_rbins(redshift_range, num_rbins, use_log10=use_log10)
        self.set_mbins(metallicity_range, num_mbins, use_log10=use_log10)

        # 3D table: axes are [temperature, redshift, metallicity]
        self.etable = np.zeros([self.num_tbins, self.num_rbins, self.num_mbins])

        # Outer loop over redshift: re-initialise the CIE session each time
        # because the redshifted energy grid passed to pyatomdb changes.
        for ir, rbin in enumerate(tqdm(self.rbins, desc="redshift")):
            z = rbin - 1.0  # convert stored (1+z) back to z
            self.setup_cie_xray_spectrum(z, nH=nH)
            for it, t in enumerate(self.tbins):
                for im, m in enumerate(self.mbins):
                    self.etable[it, ir, im] = self.compute_xray_emissivity(
                        t, m, z, nH
                    )

    # ---------------------------------------------------------------------- #
    # I/O
    # ---------------------------------------------------------------------- #

    def save_emissivity_table(self, filename, overwrite=True):
        """
        Save the 3D emissivity table to an HDF5 file.

        Parameters
        ----------
        filename : str
            Output filename. ``.hdf5`` is appended if not already present.
        overwrite : bool, optional
            Overwrite an existing file. Default True.
        """
        if self.etable is None:
            logging.error("Emissivity table not set up yet!")
            raise RuntimeError("Emissivity table not set up yet!")

        filepath = Path(filename)
        if filepath.suffix != ".hdf5":
            filepath = filepath.with_name(filepath.name + ".hdf5")
        if filepath.exists() and not overwrite:
            print("File already exists! Not saving file.")
            return

        print(f"Saving emissivity table as {filepath}")
        with h5py.File(filepath, "w") as f:
            f.create_dataset("emissivity_table", data=self.etable, dtype="f8")
            f.create_dataset("temperature_array", data=self.tbins, dtype="f8")
            f.create_dataset("redshift_array",    data=self.rbins, dtype="f8")
            f.create_dataset("metallicity_array", data=self.mbins, dtype="f8")
            f.create_dataset("energy_array",      data=self.ebins, dtype="f8")

    def read_emissivity_table(self, filename):
        """
        Read a previously saved 3D emissivity table from an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to the HDF5 file.
        """
        try:
            f = h5py.File(filename, "r")
        except OSError as err:
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), filename
            ) from err
        with f:
            self.etable = f["emissivity_table"][:]
            self.tbins  = f["temperature_array"][:]
            self.rbins  = f["redshift_array"][:]
            self.mbins  = f["metallicity_array"][:]
            self.ebins  = f["energy_array"][:]
            self.num_tbins = self.etable.shape[0]
            self.num_rbins = self.etable.shape[1]
            self.num_mbins = self.etable.shape[2]
            self.num_ebins = len(self.ebins)
            self.emin = self.ebins[0]
            self.emax = self.ebins[-1]

    # ---------------------------------------------------------------------- #
    # Interpolation
    # ---------------------------------------------------------------------- #

    def return_interpolated_emissivity(
        self,
        temperature,
        redshift,
        metallicity,
        method="linear",
        use_log10=True,
    ):
        """
        Return emissivity from the interpolated 3D table.

        Inputs may be scalars or arrays of broadcastable shapes.
        Values outside the tabulated range are clamped to the nearest edge.

        Parameters
        ----------
        temperature : float or ndarray
            Plasma temperature in keV.
        redshift : float or ndarray
            Source redshift z.
        metallicity : float or ndarray
            Metallicity in solar units.
        method : str, optional
            Interpolation method passed to scipy.interpolate.interpn.
            Default 'linear'. Use 'splinef2d' only for 2D tables; for 3D
            use 'linear' or 'nearest'.
        use_log10 : bool, optional
            Interpolate in log space (recommended). Default True.

        Returns
        -------
        emissivity : float or ndarray
            Emissivity in erg s^-1 cm^3 or ph s^-1 cm^3 depending on
            use_energy_unit. Negative values are set to zero.
        """
        # Build the grid axes in the interpolation space
        if use_log10:
            ltbins = np.log10(self.tbins)
            lrbins = np.log10(self.rbins)           # stored as (1+z)
            lmbins = np.log10(self.mbins)
        else:
            ltbins = self.tbins
            lrbins = self.rbins
            lmbins = self.mbins

        points = (ltbins, lrbins, lmbins)

        # Broadcast so that e.g. an array of temperatures can be used with a
        # scalar redshift and metallicity
        scalar_input = all(np.ndim(v) == 0 for v in (temperature, redshift, metallicity))
        temperature, redshift, metallicity = np.broadcast_arrays(
            np.asarray(temperature, dtype=float),
            np.asarray(redshift, dtype=float),
            np.asarray(metallicity, dtype=float),
        )
        shape = temperature.shape

        # Transform query coordinates (redshift bins are stored as 1+z)
        if use_log10:
            lt = np.log10(temperature).ravel()
            lr = np.log10(1.0 + redshift).ravel()
            lm = np.log10(metallicity).ravel()
        else:
            lt = temperature.ravel()
            lr = (1.0 + redshift).ravel()
            lm = metallicity.ravel()

        # Clamp to grid boundaries
        lt = np.clip(lt, ltbins[0], ltbins[-1])
        lr = np.clip(lr, lrbins[0], lrbins[-1])
        lm = np.clip(lm, lmbins[0], lmbins[-1])

        xi = np.column_stack([lt, lr, lm])
        result = interpn(
            points,
            self.etable,
            xi,
            method=method,
            bounds_error=False,
            fill_value=0.0,
        )

        result = np.clip(result, 0.0, None)

        if scalar_input:
            return float(result[0])

        return result.reshape(shape)