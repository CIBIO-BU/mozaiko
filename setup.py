from setuptools import find_packages, setup

setup(
    name="mozaiko",
    version="0.1.4",
    author="Camila Babo",
    author_email="camila.babo@cibio.up.pt",
    url="https://github.com/CIBIO-BU/mozaiko",
    packages=find_packages(),
    package_data={"": ["*"]},
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "biopython",
        "tqdm",
    ],
    entry_points={
        "console_scripts": ["mozaiko=src.mozaiko.mozaiko:main"],
    },
)