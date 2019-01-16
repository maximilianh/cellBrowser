Advanced Topics
---------------

To update only the javascript files and re-create the index.html, you can use the command line tool cbUpgrade.

To add Google Analytics tracking to your cell browser, create a file ``.cellbrowser.conf`` in your home directory
and add a line like this::

    gaTag = "UA-11231232-1"

Then cbBuild or cbUpgrade and your index.html should contain the Google Analytics tracking code.

The html directory can be defined in all tools with the option ``-o``. If that
becomes cumbersome, you can also permanently set it through the environment
variable CBOUT or by adding a line like this to ``~/.cellbrowser.conf``::

    htmlDir = "/data/www/cb/"
