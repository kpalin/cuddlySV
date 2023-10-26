# coding=utf-8

from setuptools import setup

with open("README.md") as f:
    readme = f.read()

setup(
    name="cuteSV",
    version="3.0.0",
    description="Long-read-based human genomic structural variation detection with cuteSV",
    author="Jiang Tao",
    author_email="tjiang@hit.edu.cn",
    url="https://github.com/tjiangHIT/cuteSV",
    license="MIT",
    packages=["cuteSV"],
    package_dir={"": "src/"},
    data_files=[("", ["LICENSE"])],
    scripts=["src/cuteSV/cuteSV"],
    # long_description = LONG_DESCRIPTION,
    long_description=readme,
    long_description_content_type="text/markdown",
    zip_safe=False,
    install_requires=[
        "pysam",
        "Biopython",
        "Cigar",
        "numpy",
        "pyvcf",
        'importlib-metadata >= 1.0 ; python_version < "3.8"',
    ],
)
