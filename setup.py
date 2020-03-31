import setuptools
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

setup(
    name="mdwc3",
    description="A Python library for ab-initio molecular dynamics simulations.",
    version="0.3",
    author="Arturo Hernandez",
    author_email="my@email.com",
    url="https://github.com/romerogroup/mdwc3",
    download_url="https://github.com/romerogroup/mdwc3/archive/0.3.tar.gz",
    packages=['mdwc', 'mdwc.software_tools', 'mdwc.MD_suite'],
    license="LICENSE.txt",
    scripts=['bin/mdwc3'],
    install_requires=["numpy>=1.11.1",],
    ext_modules=[
        Extension(
            name='mdwc.MD_suite.MD_suite',
            sources=['mdwc/MD_suite/MD_suite.f90'],
            extra_compile_args=["-I.", "-O3", "-sse2"],
        )
    ],
)
