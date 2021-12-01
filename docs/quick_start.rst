Quick Start Guide
----

These steps will walk through downloading and setting a minimal cell browser using an example dataset. 
It is assumed that you are carrying out these steps locally on your computer. 

It is based on `Nowakowski et al. 2017. <https://science.sciencemag.org/content/358/6368/1318.long>`_
and the corresponding https://cortex-dev.cells.ucsc.edu dataset. The expression matrix
only includes 100 genes, but it does show off many of the features of the cell browser.



Step 1: Download the example dataset
^^^^

Download the data using curl::

  curl -ks https://cells.ucsc.edu/downloads/samples/mini.tgz 

This file contains all the pieces necessary to create a simple cell browser. 

Step 2: Unpack the dataset
^^^^

Next, we need to take this file and unpack it's contents::

  tar xvz mini.tgz

This will create a new ``mini`` directory contain the cell browser files. 

Step 3: Enter the `mini` directory
^^^^

Now, enter the ``mini`` directory so that we can build the cell browser next::

  cd mini
  
Step 4: Build the cell browser
^^^^

Then, build the cell browser:: 

  cbBuild -o ~/public_html/cells/ -p 8888

The ``-o`` option specifies the output directory and ``-p`` defines the port to used
for the webserver. A bunch of information will be displayed on the screen, mostly just reporting
its progress in building the cell browser. 

Step 5: Check it out!
^^^^

Finally, point your web browser to http://localhost:8888 to view your minimal cell browser.

Step 6: Stop webserver
^^^^

When you are done playing with your minimal cell browser, stop the web server started by cbBuild using ``Ctrl-C``. 

..
  Commenting this out for now
  ----

  Building your own Cell Browser
  ^^^^

  The next page, "Setup Your Own", will describe the process of building a cell browser for your own dataset. 
