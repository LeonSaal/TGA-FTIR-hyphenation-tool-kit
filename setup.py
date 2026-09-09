import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="TGA_FTIR_tools",
    version="v2.2",
    author="Leon Saal",
    author_email="mail.leon.saal@gmail.com",
    description="A package for handling hyphenated TGA and FTIR data. Includes basic plotting as well as as advanced deconvolution of EGA data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    install_requires=[
        "pandas[excel, performance, plot, output-formatting, computation, output-formatting]>=3.0.3",
        "scikit-learn>=1.9.0",
        "requests>=2.34.2",
        "pint>=0.25.3",
        "ipykernel>=7.3.0",
        "seaborn>=0.13.2",
        "molmass>=2026.6.9",
        "pint-pandas>=0.8.0",
    ],
    python_requires=">=3.10.4",
    package_data={"settings": ["*.ini", "*.xlsx", "*.json"]},
)
