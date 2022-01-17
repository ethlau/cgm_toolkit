from setuptools import setup
from setuptools import find_packages


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
    ],
    package_dir={"":"src"},
    packages=find_packages(where="src"),
    install_requires=[
        'astropy',
        'h5py',
        'matplotlib',
        'numpy',
        'scipy',
        'pyatomdb'],
    python_requires=">=3.6"
)
