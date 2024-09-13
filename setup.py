from setuptools import setup, find_packages

setup(
    name='bioinfutils',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'bioinfutil=bioinfutils.cli:main',
        ],
    },
    install_requires=[
        Bio, argparse, pandas, os, sys
    ],
    #test_suite='tests'  # Add tests directory for automatic test discovery
)