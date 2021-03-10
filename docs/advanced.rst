Advanced Topics
---------------

Cell browser updates and web server configuration
^^^^

New features and bug fixes are being added to the UCSC Cell Browser software all the time. You can update the javascript files and re-create the index.html using the command line tool ``cbUpgrade``. You need to use the ``--code`` option and the ``-o`` to specify the output directory, e.g. 
``cbUpgrade --code -o /var/www/cellbrowser/``.

Your web server should support byte-range requests. This isn't important for smaller,
but for datasets with files larger than 30MB, a warning message
will be shown once. Most web servers and web hosters support them by default.
For Apache, byte-range requests are enabled by default but may need to be
activated in some installations of nginx.

Default output directory for building cell browsers
^^^^^

The html directory can be defined in all tools with the option ``-o``. If that
becomes cumbersome, you can also permanently set it through the environment
variable CBOUT (e.g. in your ``~/.bashrc``) or by adding a line like this to ``~/.cellbrowser.conf``::

    htmlDir = "/data/www/cb/"

Google Analytics
^^^^

To add Google Analytics tracking to your cell browser, create a file ``.cellbrowser.conf`` in your home directory
and add a line like this::

    gaTag = "UA-11231232-1"

Then run ``cbBuild`` or ``cbUpgrade`` to rebuild your index.html, after which it
should contain your Google Analytics tracking code.

Various ``cellbrowser.conf`` configurations
^^^^

For a reference of all possible ``cellbrowser.conf`` statements, see the `example conf  <https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf>`_

If you have meta fields with very long names, you can reduce the font size. Configure them like this:

    metaOpt = {'Cluster_field' : {'fontSize': '10px'}}

In your meta.tsv, you can have URLs to images. These will be shown on mouse over in the left annotation bar. 

If you set the default coloring field to 'None' (without the quotes), then there is no coloring at all when the
cell browser starts.

To change the coloring/label field automatically when the user activates some coordinates (layout), use the option
"colorOnMeta" to specify the field:: 

    coords=[
        {"file":"tsne.coords.tsv", "shortLabel":"t-SNE on WGCNA"},       
        # you can force coloring of some other meta data field when a layout is changed to another one
        {"file":"subset.coords.tsv", "shortLabel":"neural cells", colorOnMeta="neuralCluster"},
    ]

In very rare cases, it can be necessary to tell ``cbBuild`` that the numbers in the matrix are floating point numbers. 
The setting looks like this::

    matrixType = "float"
