from setuptools import setup, find_packages
import os
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8').read()
setup(
    name = "aRST-chem",
    version = "2.0.2",
    author = "Ying Chen",
    author_email = "ying.chen@pc.rwth-aachen.de",
    description = ("an Automated Reaction Search Tool"),
    license = "MIT",
    url = "https://git.rwth-aachen.de/bannwarthlab/aRST",
    packages=find_packages(),
    exclude_package_data={'': ['.gitignore']},
    include_package_data=True,
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
    ],
    install_requires=[
        "billiard~=4.2.1",
        "matplotlib~=3.7.5",
        "molbar~=1.1.3",
        "networkx~=3.1",
        "numpy~=1.24.4",
        "pandas~=2.0.3",
        "read~=0.0.2",
        "scipy~=1.10.1",
        "setuptools~=70.0.0",
        "toml~=0.10.2"
    ],
    entry_points={
        'console_scripts': [
            'aRST = aRST.script.main:main',
        ]
    }
)
