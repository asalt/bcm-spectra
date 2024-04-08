from setuptools import setup, find_packages

setup(
    name="bcm_spectra",
    version="0.1.0",
    # packages=find_packages(),
    package_dir={"": "src"},
    author="Alexander Saltzman",
    author_email="",
    description="A Python package for processing and analyzing mass spectrometry data",
    # long_description=open('README.md').read(),
    # long_description_content_type='text/markdown',
    url="https://github.com/asalt/bcm-spectra",
    license="MIT",
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "pyteomics",
        "spectrum_utils",
        "tqdm",
        #'psims'
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        # Add any other relevant classifiers
    ],
    python_requires=">=3.7",
    entry_points={"console_scripts": ["bcm_spectra=bcm_spectra.cli:main"]},
    include_package_data=True,
)
