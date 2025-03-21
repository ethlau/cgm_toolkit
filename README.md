# CGM Toolkit

[![DOI](https://zenodo.org/badge/891676081.svg)](https://doi.org/10.5281/zenodo.14714182) [![License](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Overview

The CGM Toolkit is used for analytic modeling observables of the circumgalactic medium in X-ray, Sunyaev-Zeldovich (SZ) effect, and UV absorption lines. 

## Types and Formats of Files Uploaded

- **Jupyter Notebooks:** `.ipynb` files containing code, documentation, and usage examples.
  - `example.ipynb`: Provides a usage example of the toolkit.
- **Python Scripts:** `.py` files containing the core functionality and utilities for data processing and analysis.

## Tools and Codes

- **PyAtomDB:** Used for X-ray calculations, utilizing the APEC model. More details can be found in the [PyAtomDB documentation](https://atomdb.readthedocs.io/en/master/).
- **Installation:** Install the toolkit using `pip` with the following command:
  ```bash
  pip install --user -e .
  ```

## Usage of Uploaded Scripts and Codes

- **example.ipynb:** Demonstrates how to use the toolkit for analyzing CGM data, including data preprocessing, visualization, and analysis.

## References

The CGM Toolkit has been used in the following papers:
- Lau et al 2025: [X-raying CAMELS: Constraining Baryonic Feedback in the Circum-Galactic Medium with the CAMELS simulations and eRASS X-ray Observations](https://ui.adsabs.harvard.edu/abs/2024arXiv241204559L/abstract)
- Hernández-Martínez et al 2025: [Cosmological and Astrophysical Parameter Inference from Stacked Galaxy Cluster Profiles Using CAMELS-zoomGZ](https://ui.adsabs.harvard.edu/abs/2025ApJ...981..170H/abstract)
- Singh et al 2024: [Comparison of models for the warm-hot circumgalactic medium around Milky Way-like galaxies](https://ui.adsabs.harvard.edu/abs/2024MNRAS.532.3222S/abstract)


## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes or enhancements.

## Contact

For any questions or inquiries, please contact [ethlau](https://github.com/ethlau).

Simple code to compute X-ray, SZ, and UV absorption line observables for given temperature and metallicity of the halo gas 
Uses the APEC model in [PyAtomDB](https://atomdb.readthedocs.io/en/master/) for the X-ray calculations. 
Usage example is provided in the ipython notebook `example.ipynb`. 

## Install

In the download folder, use `pip` to install

```bash
$ pip install --user -e .
```
