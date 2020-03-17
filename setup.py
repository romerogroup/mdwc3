import setuptools
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

setup(
    name="mdwc",
    description="A Python library for ab-initio molecular dynamics simulations.",
    version="0.2",
    author="Arturo Hernandez",
    author_email="my@email.com",
    packages=['mdwc', 'mdwc.software_tools', 'mdwc.MD_suite'],
    scripts=['bin/mdwc3'],
    install_requires=["numpy>=1.17.2",],
    ext_modules=[
        Extension(
            name='mdwc.MD_suite.MD_suite',
            sources=['mdwc/MD_suite/MD_suite.f90'],
            extra_compile_args=["-I.", "-O3"],
        )
    ],
)
