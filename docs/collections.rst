Dataset Collections
-------------------

The Cell Browser allows you to group related datasets into single-level collections. Collections will displayed
in your cell browser alongside all of your other datasets; when you open a collection, it will
show you all of the datasets within it.

First, to enable the collections feature, you must add a single line pointing to
the directory where all of your single-cell data lives to your ``~/.cellbrowser.conf``::

    collDir='/celldata/'

Next, for each dataset you would like to be a part of a collection,
add a ``collections`` line to the dataset's ``cellbrowser.conf``, such as::

    collections = ["organoids"]

This will create a single collection in your cell browser named ``organoids``. You can
specify multiple collection names separated by a comma (e.g. ``["organoids", "human"]``).

Then, describe the menu entry for the collection by placing a ``cellbrowser.conf`` for it somewhere within
your ``collDir``. This minimal ``cellbrowser.conf`` file only needs to contains the ``name``,
``shortLabel`` and ``tags`` settings::

   mkdir -p /celldata/organoids
   cd /celldata/organoids
   echo 'name="organoids"' > cellbrowser.conf
   echo 'shortLabel="Brain Organoids"' >> cellbrowser.conf
   echo 'tags=["10x"]' >> cellbrowser.conf

Now you can describe your collection as discussed under the **Describing
datasets** section. Put the ``desc.conf`` file into the same directory as the
``cellbrowser.conf`` you just created.

Now run ``cbBuild`` for each of the datasets that you would like to be in the collection.
If you view your cell browser on the web, you should see this new collection present.
Additionally, when viewing a dataset in a collection, you can move quickly between
it and other datasets in the same collection using the "Collection" dropdown menu.

