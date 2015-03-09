from setuptools import setup, find_packages

setup(name="betasearch",
    version="0.1",
    description="a python module for the fast indexing and querying of protein beta-sheet substructures.",
    author="Kian Ho",
    author_email="hui.kian.ho@gmail.com",
    url="http://www.github.com/kianho/betasearch",
    packages=["betapy", "ptgraph2"],
    package_dir={"" : "betasearch"},
    install_requires=["numpy", "biopython", "Whoosh"],
    setup_requires=["numpy"],
    zip_safe=False,
    keywords=["protein", "protein structure", "bioinformatics",
        "computational biology", "amino acids", "beta-sheets",
        "substructure search", "indexing", "querying", "matrices"]
)
