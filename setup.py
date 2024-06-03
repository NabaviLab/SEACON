import os
import subprocess

from setuptools import setup, find_packages
from setuptools.command.install import install

class CustomInstallCommand(install):
    """Customized setuptools install command to run initialize script after install."""
    def run(self):
        install.run(self)
        subprocess.run(['./initialize.sh'])

setup(
    name='seacon',
    version='1.1.0',
    author='Samson Weiner',
    author_email='samson.weiner@uconn.edu',
    description='SEACON (Single-cell Estimation of Allele-specific COpy Numbers) is a tool for allele-specific copy number profiling from single-cell DNA-sequencing data..',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/NabaviLab/SEACON',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'seacon=seacon.core:main',
        ],
    },
    cmdclass={
        'install': CustomInstallCommand,
    }
)