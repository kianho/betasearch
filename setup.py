from setuptools import setup, find_packages

setup(name="betasearch",
    version="0.1",
    description="a python module for the fast indexing and querying of protein beta-sheet substructures.",
    author="Kian Ho",
    author_email="hui.kian.ho@gmail.com",
    url="http://www.github.com/kianho/betasearch",
    license="MIT",
    packages=find_packages("betasearch"),
    package_dir={"" : "betasearch"},
    zip_safe=False,
    keywords=["protein", "protein structure", "bioinformatics",
        "computational biology", "amino acids", "beta-sheets",
        "substructure search", "indexing", "querying", "matrices"]
)
