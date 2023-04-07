import numpy as np
from pathlib import Path

ion_frac_file_path = Path(__file__).parent / "../data/ion_frac_CIE.txt"

class IonFrac() :

    def __init__(self, ion_table=ion_frac_file_path):

        self.ion_frac_table = np.genfromtxt(ion_table,delimiter=',', names=True)

    def return_ion_fraction (self, temperature, ion) :

        table = self.ion_frac_table

        input_temperature = table['T']

        input_ion_frac = table[ion]

        frac = 10**np.interp(np.log10(temperature), np.log10(input_temperature), np.log10(input_ion_frac))

        return frac

