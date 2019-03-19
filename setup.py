import setuptools
import versioneer

with open("pypi/README.md", "r") as fh:
    long_description = fh.read()

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
        'cbImportSeurat2 = cellbrowser.seurat:cbImportSeurat2Cli',
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

    zip_safe = False, # do not allow install as an egg, we need to access the .html, .js, .css, etc files
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: JavaScript"
    ],
)
