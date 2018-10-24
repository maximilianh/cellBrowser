import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cellbrowser",
    version="0.4.6",
    python_requires='>=2.5',
    author="Maximilian Haeussler",
    author_email="max@soe.ucsc.edu",
    url="https://github.com/maximilianh/cellBrowser",
    description="UCSC Cellbrowser for single cell data. Includes command line builder and dataset converter tools for Seurat, Scanpy and Cellranger.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages("../src/cbPyLib/"),
    package_dir={'':'../src/cbPyLib/'},   # tell distutils packages are under src
    include_package_data=True,  # use MANIFEST.in for non-python files
    #package_data={
        #'cellbrowser': ['cbWeb/js/*.js', 'cbWeb/html/*', 'cbWeb/ext/*']
    #},
    entry_points={
    'console_scripts': [
        'cbBuild = cellbrowser.cellbrowser:cbBuildCli',
        'cbScanpy = cellbrowser.cellbrowser:cbScanpyCli',
        'cbTool = cellbrowser.convert:cbToolCli',
        'cbUpgrade = cellbrowser.cellbrowser:cbMake_cli',
        'cbGuessGencode = cellbrowser.guessgenes:cbGuessGencodeCli',
        'cbMarkerAnnotate = cellbrowser.geneinfo:cbMarkerAnnotateCli'
    ]
    },
    #package_data={
        # If any package contains *.txt files, include them:
        #'': ['*.txt'],
        # And include any *.dat files found in the 'data' subdirectory
        # of the 'mypkg' package, also:
        #'mypkg': ['data/*.dat'],
    #}

    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: JavaScript"
    ],
)
