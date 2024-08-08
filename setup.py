from setuptools import setup, find_packages

setup(
    name='dnaquaimg',
    version='0.1.0',
    author='Camila Babo',
    author_email='camila.babo@cibio.up.pt',
    url='https://github.com/CIBIO-BU/DNAquaIMG',
    packages=find_packages(),
    install_requires=[
        'numpy==2.0.0',
        'pandas==2.2.2'
    ],
)