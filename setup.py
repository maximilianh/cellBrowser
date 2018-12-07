import setuptools

with open("pypi/README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cellbrowser",
    version="0.4.24",
    license="GPL 3",
    python_requires='>=2.5',
    author="Maximilian Haeussler",
    author_email="max@soe.ucsc.edu",
    url="https://github.com/maximilianh/cellBrowser",
    description="UCSC Cellbrowser for single cell data. Includes tab-sep command line importer, dataset converter tools for Seurat, Scanpy and Cellranger and two basic pipelines to process an expression matrix.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages("src/cbPyLib/"),
    package_dir={'':'src/cbPyLib/'},   # tell distutils packages are under src
    include_package_data=True,  # use MANIFEST.in for non-python files
    #package_data={
        #'cellbrowser': ['cbWeb/js/*.js', 'cbWeb/html/*.html', 'cbWeb/ext/*']
    #},
    #scripts=['src/cbScanpy'], # want to force python3 as the executable for cbScanpy
    entry_points={
    'console_scripts': [
        'cbBuild = cellbrowser.cellbrowser:cbBuildCli',
        'cbScanpy = cellbrowser.cellbrowser:cbScanpyCli',
        'cbSeurat = cellbrowser.seurat:cbSeuratCli',
        'cbTool = cellbrowser.convert:cbToolCli',
        'cbUpgrade = cellbrowser.cellbrowser:cbMake_cli',
        'cbGuessGencode = cellbrowser.guessgenes:cbGuessGencodeCli',
        'cbMarkerAnnotate = cellbrowser.geneinfo:cbMarkerAnnotateCli',
        'cbImportScanpy = cellbrowser.convert:cbImportScanpyCli',
        'cbImportCellranger = cellbrowser.convert:cbCellrangerCli'
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
