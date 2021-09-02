#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = ['click',
                'decorator',
                'flask',
                'ftputil',
                'lxml',
                'matplotlib',
                'networkx',
                'pytest',
                'pandas',
                'requests',
                'scipy',
                'tox',
                'tqdm']

test_requirements = ['pytest>=3', ]

setup(
    author="Lauren D., Rebeca F., Simon M., Maren P.",
    author_email='laurendelong21@gmail.com, maren.philipps@uni-bonn.de',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Cross Verification of Protein Information",
    entry_points={
        'console_scripts': [
            'dbi=dbinspector.cli:cli',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='dbinspector',
    name='dbinspector',
    packages=find_packages(include=['dbinspector', 'dbinspector.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
