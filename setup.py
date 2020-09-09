from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='tsBNgen',
    version='1.0.0',
    author='Manie Tadayon',
    author_email='manitadayon@ucla.edu',
    description='Generate time series data from an arbitrary Bayesian network',
    packages=['tsBNgen'],
    license='MIT',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/manitadayon/tsBNgen',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)