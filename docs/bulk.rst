Bulk downloads
--------------

The dataset hierarchy
^^^^^^^^^^^^^^^^^^^^^

There is no search API yet but a systematic way to download::
1) http://cells.ucsc.edu/dataset.json contains the list of top-level datasets in the "datasets" attribute.
2) every object in the "datasets" list is a child dataset, note the "name" attribute.
3) every dataset has a <datasetname>/dataset.json file again, e.g. http://cells.ucsc.edu/adultPancreas/dataset.json
4) datasets with a "datasets" attribute have children themselves, e.g. http://cells.ucsc.edu/xena/dataset.json

Rsync download
^^^^^^^^^^^^^^

If you want to download all or selected datasets, use our ``rsync`` server. It is still being tested,
please email us for more info.
