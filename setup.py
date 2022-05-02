from setuptools import setup, find_packages
from shutil import which


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


if not is_tool("seqkit"):
    print("seqkit must be installed for this package to work!")
    exit()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="barcode_demultiplex",
    version="0.0.1",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    packages=["barcode_demultiplex"],
    py_modules=[
        "barcode_demultiplex/run",
    ],
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "barcode_demultiplex = barcode_demultiplex.run:cli"
        ]
    },
)
