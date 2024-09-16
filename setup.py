from setuptools import setup, find_packages

setup(
    name='bioinfutils',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'bioinfutils=bioinfutils.cli:main',
        ],
    },
    install_requires=[
        'Bio', 'argparse', 'pandas'
    ]
    #test_suite='tests'  # Add tests directory for automatic test discovery
)