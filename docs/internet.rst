Putting it onto the internet
----------------------------

To deploy the result of ``cbBuild`` onto a webserver, simply copy all files and
directories in the html output directory (in the examples on this site, that's
usually *~/public_html/cells*) to an empty directory on a webserver and point
your web browser to it. Many universities give their members
webspace, sometimes in a directory called ~/public_html or on a special server.
If you do not have one, contact us or use Cyverse, Amazon S3, Goole Cloud
Storage, Microsoft Azure etc. to host your files. You cannot use online backup
solutions like Dropbox, Box.com, iCloud OneDrive or Google Drive, they
intentionally are not webservers. You can always send the result to cells@ucsc.edu, 
we are happy to add it to our cell browser site at cells.ucsc.edu.

To add more datasets, go to the other data directories and run cbBuild
there, with the same output directory. cbBuild will then modify the index.html
in the output directory to show all datasets. Note that the directory that you
provide via -o (or the CBOUT environment variable) is the html directory. The
data for each individual dataset will be copied into subdirectories under this
html directory, one directory per dataset.

Instead of specifying "-o" all the time, you can also add a line like this to
your ~/.bashrc to point to your html directory::
 
    export CBOUT=/var/www

The -p 8888 is optional. A more permanent alternative to the -p option is to
run a webserver on your machine and build directly into its web directory.

On a Mac you can use the Apache that ships with OSX::

    sudo /usr/sbin/apachectl start
    sudo cbBuild -o /Library/WebServer/Documents/cells/

Then you should be able to access your viewer at http://localhost/cells

On Linux, you would install Apache2 (with 'sudo yum install htppd' or 'sudo apt-get install
apache2') and use the directory /var/www/ instead::

    sudo cbBuild -o /var/www/

We hope you do not use this software on Windows. Email cells@ucsc.edu if you have to.


