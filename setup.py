from setuptools import setup, find_packages

version = {}
with open("dnachisel/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="dnachisel",
    version=version["__version__"],
    author="Zulko",
    description="Optimize DNA sequences under constraints.",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/DnaChisel",
    keywords=("DNA gene design codon optimization constraints " "synthetic biology"),
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    scripts=["scripts/dnachisel"],
    install_requires=[
        "numpy",
        "biopython",
        "proglog",
        "docopt",
        "flametree",
        "python_codon_tables",
    ],
    extras_require={
        "reports": [
            "pdf_reports",
            "sequenticon",
            "matplotlib",
            "dna_features_viewer",
            "pandas",
        ],
        "tests": [
            "pytest",
            "pytest-cov",
            "coveralls",
            "geneblocks",
            "genome_collector",
            "matplotlib",
            "primer3-py",
        ],
    },
)
