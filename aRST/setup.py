from setuptools import setup, find_packages
# import os
# def read(fname):
#     return open(os.path.join(os.path.dirname(__file__), fname)).read()
setup(
    name = "aRST",
    version = "1.0.0",
    author = "Ying Chen",
    author_email = "ying.chen@pc.rwth-aachen.de",
    description = ("an automated Reaction Search Tool"),
    license = "unkown",
    url = "https://git.rwth-aachen.de/bannwarthlab/reaction_discovery_python_script",
    packages=find_packages(include=['aRST']),
    exclude_package_date={'':['.gitignore']},
    # long_description=read('README'),
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
    ],
    entry_points={
        'console_scripts': [
            'aRST = aRST.script.main:main',
        ]
    }
)
