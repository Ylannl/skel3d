from setuptools import setup, find_packages
setup(
    name = "skel3d",
    version = "0.1",
    packages = find_packages(),
    install_requires = ["numpy>=1.9.1", "povi", "click", "scipy", "sklearn", "python-igraph"],
    include_package_data = True,
    # extras_require = {
    #     'LAS':  ["laspy"]
    #},
    scripts = ['scripts/segment.py']
)