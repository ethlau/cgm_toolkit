"""
Generate a 2D X-ray emissivity table at a fixed redshift z = 0.

The table holds the CIE emissivity Lambda(T, Z) in the 0.5-2 keV band (by
default), tabulated over temperature (keV) and metallicity (solar units).
It is written in the same 2D layout as the older tables in data/
(e.g. etable_05_2keV.hdf5):

    emissivity_table   shape (num_tbins, num_metal_bins)
    temperature_array  keV
    metallicity_array  Zsun
    energy_array       keV

The defaults match etable_05_2keV.hdf5:

    band         0.5 - 2.0 keV, 101 energy bins
    temperature  0.000862 - 30 keV, 101 log-spaced bins
    metallicity  0.001 - 10 Zsun, 21 log-spaced bins
    units        erg s^-1 cm^3, no Galactic absorption, no instrument response

Note: this 2D layout has no redshift axis, so it cannot be read by
XrayEmissivity.read_emissivity_table or used by cgm_profiles, which expect
3D tables (see make_etable_3D.py). Read it directly with h5py.

Requires pyatomdb with the AtomDB data files installed.

Examples
--------
    python scripts/make_etable_2D.py

    python scripts/make_etable_2D.py -o data/etable_2D_cnts.hdf5 \\
        --photons --nH 2e20 \\
        --rmf responses/rmf01_sdtq.fits --arf responses/arf01_200nmAl_sdtq.fits
"""

import argparse
from pathlib import Path

import h5py

from cgm_toolkit.xray_emissivity import XrayEmissivity

default_output = Path(__file__).parent / "../data/etable_05_2keV_2D.hdf5"


def parse_args():
    p = argparse.ArgumentParser(
        description="Tabulate a 2D (temperature, metallicity) X-ray "
        "emissivity table at z = 0.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-o", "--output", default=str(default_output),
                   help="Output HDF5 file")
    p.add_argument("--emin", type=float, default=0.5,
                   help="Lower edge of the energy band in keV")
    p.add_argument("--emax", type=float, default=2.0,
                   help="Upper edge of the energy band in keV")
    p.add_argument("--num-ebins", type=int, default=101,
                   help="Number of energy bins")
    p.add_argument("--tmin", type=float, default=0.000862,
                   help="Minimum T in keV")
    p.add_argument("--tmax", type=float, default=30.0, help="Maximum T in keV")
    p.add_argument("--num-tbins", type=int, default=101,
                   help="Number of temperature bins")
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

    output = Path(args.output)
    if output.suffix != ".hdf5":
        output = output.with_name(output.name + ".hdf5")
    if output.exists() and args.no_overwrite:
        print(f"{output} already exists! Not saving file.")
        return

    # A single redshift bin at z = 0 gives a [T, 1, Z] table.
    xray = XrayEmissivity()
    xray.tabulate_xray_emissivity(
        energy_range=[args.emin, args.emax],
        num_ebins=args.num_ebins,
        rmf=args.rmf,
        arf=args.arf,
        use_energy_unit=not args.photons,
        temperature_range=[args.tmin, args.tmax],
        redshift_range=[0.0, 0.0],
        metallicity_range=[args.metal_min, args.metal_max],
        nH=args.nH,
        num_tbins=args.num_tbins,
        num_rbins=1,
        num_mbins=args.num_metal_bins,
        use_log10=True,
    )
    etable_2d = xray.etable[:, 0, :]

    print(f"Saving 2D emissivity table as {output}")
    with h5py.File(output, "w") as f:
        f.create_dataset("emissivity_table", data=etable_2d, dtype="f8")
        f.create_dataset("temperature_array", data=xray.tbins, dtype="f8")
        f.create_dataset("metallicity_array", data=xray.mbins, dtype="f8")
        f.create_dataset("energy_array", data=xray.ebins, dtype="f8")


if __name__ == "__main__":
    main()
