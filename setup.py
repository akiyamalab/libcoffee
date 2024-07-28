from setuptools import setup

with open("README.md") as f:
    long_description = f.read()

setup(
    name="libcoffee",
    version="0.2.0",
    description="A library for compound filtering via fragment-based efficient evaluation",
    author="Keisuke Yanagisawa",
    author_email="yanagisawa@c.titech.ac.jp",
    license="MIT",
    url="https://github.com/akiyamalab/libcoffee",
    install_requires=["openbabel-wheel", "rdkit", "rdkit-stubs", "numpy", "pytest"],
    extras_require={},
    entry_points={},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    packages=["libcoffee"],
    package_data={},
    long_description=long_description,
    long_description_content_type="text/markdown",
)
