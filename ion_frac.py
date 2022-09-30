import numpy as np

class IonFrac() :

    def __init__(self, ion_table='ion_frac_CIE.txt'):

        self.ion_frac_table = np.genfromtxt(ion_table,delimiter=',', names=True)

    def return_ion_fraction (self, temperature, ion) :

        table = self.ion_frac_table

        input_temperature = table['T']


        input_ion_frac = table[ion]

        frac = np.interp(temperature, input_temperature, input_ion_frac)

        return frac

