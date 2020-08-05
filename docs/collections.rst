Dataset Hierarchies
-------------------

The Cell Browser allows you to group related datasets into a hierarchy, where
datasets are grouped into collections, like files are grouped into directories. 
When you open a collection, it will show you all of the datasets within it.

This requires your datasets to be arranged in directories on disk. Let's say
you have two directories with data files, one in directory ``data1`` and one in
directory ``data2``, each with their own cellbrowser.conf files, then these
two directories must be both subdirectories of a parent directoryn
named e.g. ``dataParent``. The names of all the datasets
are the names of their directories, not the names 
specified via the ``name`` statement in their ``cellbrowser.conf`` files anymore.
Specified names from cellbrowser.conf are ignored when dataset hierarchies are used
and are replaced with their directory names.

To enable dataset hierarchies, you only have to add a single line pointing to
the top-level parent directory where all of your single-cell data lives. 
Add a statement like the following to your ``~/.cellbrowser.conf``::

    dataRoot='/celldata/'

Alternatively, ``dataRoot`` can be set using the ``CBDATAROOT`` environment variable::

    export CBDATAROOT='/celldata/'

Then, create a "stub" cellbrowser.conf into this directory, it should only contain
a single line like ``shortLabel="some description"``. 
You can describe your collection as discussed under the **Describing
datasets** section. Put the ``desc.conf`` file into the same directory as the
``cellbrowser.conf`` you just created.
Define at least the statements ``title`` and ``description``.  They will be
shown at the top of your dataset list. This directory can be called the
top-level collection.

Arrange your dataset directories under this directory. You can add empty directories,
which will become collections, by creating dataset directories in them and put a
``cellbrowser.conf`` and ``desc.conf`` into it, e.g. like this::

   mkdir -p /celldata/organoids
   cd /celldata/organoids
   echo 'name="organoids"' > cellbrowser.conf
   echo 'shortLabel="Brain Organoids"' >> cellbrowser.conf
   echo 'tags=["10x"]' >> cellbrowser.conf

Now you can run ``cbBuild`` in each subdirectory of a collection.
in the collection.  Or you can rebuild in all subdirectories using ``cbBuild
-r``.

If you view the cell browser now using a web browser, you should see this new
collection present. When viewing a dataset in a collection, you
can move quickly to any other dataset in the same collection using the
"Collection" dropdown menu in the toolbar.

