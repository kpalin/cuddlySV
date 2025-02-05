# coding=utf-8

from setuptools import setup

with open("README.md") as f:
    readme = f.read()

setup(
    name="cuddlySV",
    version="3.0.0.dev1",
    description="Long-read-based human genomic (somatic) structural variation detection",
    author="Jiang Tao, Kimmo Palin",
    url="https://github.com/kpalin/cuddlySV",
    license="MIT",
    packages=["cuddlySV"],
    package_dir={"": "src/"},
    package_data={"": ["LICENSE"]},
    scripts=[
        # "src/cuddlySV/cuddlySV",
        "src/somatic/merge_work_dirs.sh",
        "src/somatic/add_mapping_tags.py",
        "src/somatic/make_panel_of_normals_cleaner.sh",
    ],
    entry_points={"console_scripts": ["cuddlySV=cuddlySV:cuddlySV.run"]},
    # long_description = LONG_DESCRIPTION,
    long_description=readme,
    long_description_content_type="text/markdown",
    zip_safe=False,
    python_requires=">=3",
    data_files=[
        (
            "share/cuddlySV",
            ["src/somatic/somaticSV.smk", "src/somatic/cuddly_somatic.json"],
        )
    ],
    install_requires=[
        "pysam",
        "ncls",
        "Cigar",
        "numpy",
        "pyfastx",
        "scipy",
        'importlib-metadata >= 1.0 ; python_version < "3.8"',
    ],
)
