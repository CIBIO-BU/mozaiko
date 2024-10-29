from setuptools import find_packages, setup

setup(
    name="dnaquaimg",
    version="0.1.0",
    author="Camila Babo",
    author_email="camila.babo@cibio.up.pt",
    url="https://github.com/CIBIO-BU/DNAquaIMG",
    packages=find_packages(),
    package_data={"": ["*"]},
    include_package_data=True,
    install_requires=[
        "numpy==2.0.0",
        "pandas==2.2.2",
        "matplotlib==3.9.1.post1",
        "biopython==1.78",
        "tqdm",
    ],
    entry_points={
        "console_scripts": ["mozaiko=src.mozaiko:main"],
    },
)
