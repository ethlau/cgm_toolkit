"""
Generate the 3D X-ray emissivity table used by cgm_profiles.

The table holds the CIE emissivity Lambda(T, z, Z) in an observer-frame
energy band, tabulated over temperature (keV), redshift and metallicity
(solar units). The defaults reproduce data/etable_05_2keV_3D.hdf5:

    band         0.5 - 2.0 keV, 101 energy bins
    temperature  0.001 - 20 keV, 101 log-spaced bins
    redshift     0 - 4, 51 bins log-spaced in (1 + z)
    metallicity  0.001 - 10 Zsun, 21 log-spaced bins
    units        erg s^-1 cm^3, no Galactic absorption, no instrument response

Requires pyatomdb with the AtomDB data files installed. The full default
table takes a while to compute (51 x 101 x 21 spectra).

Examples
--------
Rebuild the default table:

    python scripts/make_etable_3D.py

Photon counts through an instrument response, with Galactic absorption:

    python scripts/make_etable_3D.py -o data/etable_cnts_3D.hdf5 \\
        --photons --nH 2e20 \\
        --rmf responses/rmf01_sdtq.fits --arf responses/arf01_200nmAl_sdtq.fits
"""

import argparse
from pathlib import Path

from cgm_toolkit.xray_emissivity import XrayEmissivity

default_output = Path(__file__).parent / "../data/etable_05_2keV_3D.hdf5"


def parse_args():
    p = argparse.ArgumentParser(
        description="Tabulate a 3D (temperature, redshift, metallicity) "
        "X-ray emissivity table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-o", "--output", default=str(default_output),
                   help="Output HDF5 file")
    p.add_argument("--emin", type=float, default=0.5,
                   help="Lower edge of the observer-frame band in keV")
    p.add_argument("--emax", type=float, default=2.0,
                   help="Upper edge of the observer-frame band in keV")
    p.add_argument("--num-ebins", type=int, default=101,
                   help="Number of energy bins")
    p.add_argument("--tmin", type=float, default=0.001, help="Minimum T in keV")
    p.add_argument("--tmax", type=float, default=20.0, help="Maximum T in keV")
    p.add_argument("--num-tbins", type=int, default=101,
                   help="Number of temperature bins")
    p.add_argument("--zmin", type=float, default=0.0, help="Minimum redshift")
    p.add_argument("--zmax", type=float, default=4.0, help="Maximum redshift")
    p.add_argument("--num-zbins", type=int, default=51,
                   help="Number of redshift bins")
    p.add_argument("--metal-min", type=float, default=0.001,
                   help="Minimum metallicity in Zsun")
    p.add_argument("--metal-max", type=float, default=10.0,
                   help="Maximum metallicity in Zsun")
    p.add_argument("--num-metal-bins", type=int, default=21,
                   help="Number of metallicity bins")
    p.add_argument("--nH", type=float, default=0.0,
                   help="Galactic hydrogen column density in cm^-2 (tbabs)")
    p.add_argument("--rmf", default=None, help="RMF response file")
    p.add_argument("--arf", default=None, help="ARF response file")
    p.add_argument("--photons", action="store_true",
                   help="Tabulate photon counts instead of energy (erg)")
    p.add_argument("--no-overwrite", action="store_true",
                   help="Do not overwrite an existing output file")
    return p.parse_args()


def main():
    args = parse_args()

    xray = XrayEmissivity()
    xray.tabulate_xray_emissivity(
        energy_range=[args.emin, args.emax],
        num_ebins=args.num_ebins,
        rmf=args.rmf,
        arf=args.arf,
        use_energy_unit=not args.photons,
        temperature_range=[args.tmin, args.tmax],
        redshift_range=[args.zmin, args.zmax],
        metallicity_range=[args.metal_min, args.metal_max],
        nH=args.nH,
        num_tbins=args.num_tbins,
        num_rbins=args.num_zbins,
        num_mbins=args.num_metal_bins,
        use_log10=True,
    )
    xray.save_emissivity_table(args.output, overwrite=not args.no_overwrite)


if __name__ == "__main__":
    main()
