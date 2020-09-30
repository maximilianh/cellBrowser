Putting it onto the internet
----------------------------

Basics
^^^^

Deploying your cell browser on the web is as simple as copying the output of ``cbBuild``,
including all files and directories, into to an empty directory on a web 
server. The cell browser should be able to be deployed on almost any web server, including:

* One provided by your university (often through a public_html directory in your home directory)
* Github.io provides free hosting for 1 GB of webspace and it is fast enough, see https://sansomlab.github.io/
* You can also use commercial cloud providers, such as Amazon S3, Google Cloud Storage, Microsoft Azure, though they will need a credit card.

Please consider also sending the output files to cells@ucsc.edu, 
we are more than happy to add it to our public `Cell Browser <cells.ucsc.edu>`_ website.
You can choose a dataset prefix, and then can add myDataset.cells.ucsc.edu to your manuscript.
Unfortunately, online backup solutions such as Dropbox, Box.com, iCloud, OneDrive or Google
Drive will not work; they are intentionally designed to not be usable as web servers.

Adding multiple datasets to your cell browser
^^^^^

To add more datasets to the same cell browser, navigate to the other data directories and run cbBuild
there with the same output directory. cbBuild will then modify the index.html
in the output directory to show all datasets. Note that the directory that you
provide via -o (or the CBOUT environment variable) is the html directory. The
data for each individual dataset will be copied into subdirectories under this
html directory, one directory per dataset.

Specifying a default output directory for ``cbBuild``
^^^^^

The output directory for ``cbBuild`` can be controlled using environment or .conf variables. 
This allows you to run ``cbBuild`` in a directory without needing to specify an output
directory using the "-o" option.

To control the output directory using an environment variable, add the following line to
your ~/.bashrc to point to your html directory::
 
    export CBOUT=/var/www

Replace ``/var/www`` with whatever you want the default output directory to be.

Alternatively, you can create a file called ``.cellbrowser.conf`` in your home directory
and assign a value to htmlDir::

    echo 'htmlDir = "/var/www"' >> ~/.cellbrowser.conf


Again, replace ``/var/www`` with your own dorectory. 

Notes on setting up a permanent cell browser on a local machine
^^^^^^

The port option, e.g. ``-p 8888``, is optional. When this option is specified,
it will start up its own web server. If you are running this on your local machine,
a more permanent alternative to the -p option is to run a web server on your machine
and then build directly into its web directory.

On a Mac, you can use the Apache that ships with OSX::

    sudo /usr/sbin/apachectl start
    sudo cbBuild -o /Library/WebServer/Documents/cells/

You should be able to access your viewer at http://localhost/cells

On Linux, you will need to install Apache2 (with ``sudo yum install httpd``
or ``sudo apt-get install apache2``) and use the directory ``/var/www/`` instead::

    sudo cbBuild -o /var/www/

Windows is usually not used as a web server, just use the built-in web-server via ``-p (portNumber)`` on Windows.
