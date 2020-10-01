import setuptools
import versioneer

#with open("pypi/README.md", "r") as fh:
    #long_description = fh.read()

setuptools.setup(
    name="cellbrowser",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="GPL 3",
    python_requires='>=2.5',
    author="Maximilian Haeussler",
    author_email="max@soe.ucsc.edu",
    url="https://github.com/maximilianh/cellBrowser",
    description="UCSC Cellbrowser, an interactive browser for single cell data. Includes converters and basic pipelines for text files, Seurat, Scanpy and Cellranger.",
    long_description="""The UCSC Cell Browser is an interactive browser for
    single cell data, like mRNA or ATAC-seq data. You can display
    dimensionality reductions, navigate them with the mouse or the cursor keys,
    select cells, color by genes or meta annotations and make many other
    changes. The main site runs at https://cells.ucsc.edu, but using this
    package you can also convert data yourself and build a Cell Browser HTML
    directory that can be served through any University webserver. You can try
    the Cell Browser at https://cells.ucsc.edu or read about how to convert
    data with this package on https://cellbrowser.rtfd.org.""",
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages("src/cbPyLib"),
    package_dir={'':'src/cbPyLib'},   # tell distutils packages are under src
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
        'cbUpgrade = cellbrowser.cellbrowser:cbUpgradeCli',
        'cbGuessGencode = cellbrowser.guessgenes:cbGuessGencodeCli',
        'cbMarkerAnnotate = cellbrowser.geneinfo:cbMarkerAnnotateCli',
        'cbImportScanpy = cellbrowser.convert:cbImportScanpyCli',
        'cbImportSeurat = cellbrowser.seurat:cbImportSeuratCli',
        'cbImportCellranger = cellbrowser.convert:cbCellrangerCli',
        'cbHub = cellbrowser.hubmaker:cbHubCli'
    ]
    },
    #package_data={
        # If any package contains *.txt files, include them:
        #'': ['*.txt'],
        # And include any *.dat files found in the 'data' subdirectory
        # of the 'mypkg' package, also:
        #'mypkg': ['data/*.dat'],
    #}

    zip_safe = False, # do not allow install as an egg, we need to access the .html, .js, .css, etc files
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: JavaScript"
    ],
)
