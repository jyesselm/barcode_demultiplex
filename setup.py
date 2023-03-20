from setuptools import setup, find_packages


with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="barcode_demultiplex",
    version="0.0.2",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    packages=["barcode_demultiplex"],
    py_modules=[
        "barcode_demultiplex/cli",
        "barcode_demultiplex/demuliplex",
        "barcode_demultiplex/external_cmd",
        "barcode_demultiplex/logger",
    ],
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        "console_scripts": ["barcode_demultiplex = barcode_demultiplex.cli:cli"]
    },
)
