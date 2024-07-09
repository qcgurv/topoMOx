from setuptools import setup, find_packages

setup(
    name='topomox',
    version='1.0.0',
    author='Albert Masip-Sánchez',
    author_email='albert.masip@urv.cat',
    description='Molecular topology analysis for metal oxides',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'topomox=topomox.main:main',
        ],
    },
    install_requires=[
        'numpy',
        'pandas',
        'itertools',
        'os',
        're'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: GNU AFFERO GENERAL PUBLIC LICENSE',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
