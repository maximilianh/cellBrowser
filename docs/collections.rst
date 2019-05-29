Dataset Collections
-------------------

At the moment, the Cell Browser allows you to group related datasets into a single level called 'collections'.

To group datasets into collections, you have to enable it via a statement in
``~/.cellbrowser.conf``. Add a line that points to a directory where you keep your single cell datasets e.g.::

    collDir='/celldata/'

Then, to add any dataset to one or more collections, add a line like this to its cellbrowser.conf file::

    collections = ["organoids"]

This creates a collection called ``organoids``. You now only have to define the
menu entry of this collection and describe the content of this collection. The
menu entry is defined through a minimal cellbrowser.conf file, with just name,
shortLabel and tags::

   mkdir -p /celldata/organoids
   cd /celldata/organoids
   echo 'name="organoids"' > cellbrowser.conf
   echo 'shortLabel="Brain Organoids"' >> cellbrowser.conf
   echo 'tags=["10x"]' >> cellbrowser.conf

Now you can describe your collection as explained previously under `Describing
datasets`_. Put the desc.conf file into the same directory as the
cellbrowser.conf you just created.

Now re-build the dataset that you just put into the collection with
``cbBuild``. You should see the new collection in the user interface. For every
new dataset that you want to add to the collection, run ``cbBuild`` again and
it should appear in the menu. You can move quickly between datasets of the same
collection with the "Collection" dropdown menu.

