from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cgm_toolkit",
    version="0.0.1",
    author="Erwin Lau",
    author_email="ethlau@gmail.com",
    description="Package for modeling CGM observables ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ethlau/cgm_toolkit",
    project_urls={
        "Bug Tracker": "https://github.com/ethlau/cgm_toolkit/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    packages=find_packages(),
    install_requires=[
        'astropy',
        'h5py',
        'tqdm',
        'numpy',
        'scipy',
        'pyatomdb'],
    python_requires=">=3.9",
    package_data={"cgm_toolkit": ("data/*.hdf5",)}
)
