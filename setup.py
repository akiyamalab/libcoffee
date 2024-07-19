from setuptools import setup

setup(
    name="libcoffee",
    version="0.1",
    description="A library for compound filtering via fragment-based efficient evaluation",
    author="Keisuke Yanagisawa",
    author_email="yanagisawa@c.titech.ac.jp",
    license="MIT license",
    url="https://github.com/akiyamalab/libcoffee",
    install_requires=["openbabel-wheel", "rdkit", "rdkit-stubs"],
    extras_require={
    },
    entry_points={
    },
    packages=["libcoffee"],
    package_data={},
)
