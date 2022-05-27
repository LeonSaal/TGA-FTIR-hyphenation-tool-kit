import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="TGA-FTIR-hyphenation-tool-kit",  # Replace with your own username
    version="v2.2",
    author="Leon Saal",
    author_email="mail.leon.saal@gmail.com",
    description="A package for handling hyphenated TGA and FTIR data. Includes basic plotting as well as as advanced deconvolution of FTIR data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BAMresearch/TGA-FTIR-hyphenation-tool-kit",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
