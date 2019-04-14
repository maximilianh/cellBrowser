Advanced Topics
---------------

For a reference of all cellbrowser.conf statements, see the example file https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/sampleConfig/cellbrowser.conf

To update only the javascript files and re-create the index.html, you can use the command line tool cbUpgrade.

To add Google Analytics tracking to your cell browser, create a file ``.cellbrowser.conf`` in your home directory
and add a line like this::

    gaTag = "UA-11231232-1"

Then cbBuild or cbUpgrade and your index.html should contain the Google Analytics tracking code.

The html directory can be defined in all tools with the option ``-o``. If that
becomes cumbersome, you can also permanently set it through the environment
variable CBOUT or by adding a line like this to ``~/.cellbrowser.conf``::

    htmlDir = "/data/www/cb/"

Your webserver should support byte-range requests. Smaller datasets work
without that, but for datasets with files larger than 30MB, a warning message
will be shown once. Byte-ranges are active by default in Apache but may need to
be activated in nginx. 

If you have meta fields with very long names, you can reduce the font size. Configure them like this:

    metaOpt = {'Cluster_field' : {'fontSize': '10px'}}

In your meta.tsv, you can have URLs to images. These will be shown on mouse over in the left annotation bar. 

If you set the default coloring field to 'None' (without the quotes), then there is no coloring at all when the
cell browser starts.
