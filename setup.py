import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="TGA_FTIR_tools", 
    version="v2.2",
    author="Leon Saal",
    author_email="mail.leon.saal@gmail.com",
    description="A package for handling hyphenated TGA and FTIR data. Includes basic plotting as well as as advanced deconvolution of FTIR data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BAMresearch/TGA-FTIR-hyphenation-tool-kit",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    install_requires=["scipy", "pandas", "matplotlib", "sklearn", "openpyxl", "requests"],
    python_requires="=3.10.4",
)
