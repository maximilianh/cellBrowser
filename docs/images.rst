Microscopy images
-----------------

The Cell Browser has some basic support for showing microscopy images for your
dataset. The basic model is the one of sets of images (e.g. stainings for 
different markers of the same section) arranged into categories (e.g. stages
or tissues). Some images that can be offered for download-only (e.g. TIFFS) and 
some images can be shown on the screen, these are typically reduced size
versions in JPEG format with sizes of around 1-4k pixels. You provide 
a JSON file with the description of the images. The Cell Browser then
shows an overview of the images on the "Info & Downloads" dialog box.

The default filename is ``images.json`` and has to be placed in the same
directory as the desc.conf file. Here is an example images.json file.
Contact us if you have questions or feedback on this basic system.

In the example below, there is only one image category called "main category", 
(there could be several, by duplicating the top-most object and changing
the categoryLabel in the copy). This category has two categoryImageSets, the first has the
label "DAPI SOX2 TBR2 DCX CTIP2".  This imageSet offers two files for download
and one file for display. The next imageSet has the label "AnotherImageSet" and 
also two files for download and one for display.

    [
    {
        "categoryLabel": "main category",
        "categoryImageSets": [
            {
                "setLabel": "DAPI SOX2 TBR2 DCX CTIP2",
                "downloads": [
                    {
                        "label": "all merged",
                        "file": "images/cs13-dapi-488-dcx-546-sox2-594-tbr2-647-ctip.tif"
                    },
                    {
                        "label": "DAPI",
                        "file": "images/cs13-DAPI-dcx-sox2-tbr2-ctip2.tif"
                    }
                ],
                "images": [
                    {
                        "label": "all merged",
                        "file": "images/cs13-dapi-488-dcx-546-sox2-594-tbr2-647-ctip.5000.jpg"
                        "thumb": "images/cs13-dapi-488-dcx-546-sox2-594-tbr2-647-ctip.500.jpg"
                    }
                ],
            },
            {
                "setLabel": "AnotherImageSet",
                "downloads": [
                    {
                        "label": "dicom file1",
                        "file": "images/test1.dicom"
                    },
                    {
                        "label": "dicom file2",
                        "file": "images/test2.dicom"
                    }
                ],
                "images": [
                    {
                        "label": "overview",
                        "file": "images/testMerge.jpg"
                        "thumb": "images/testMerge.thumb.jpg"
                    }
                ]
            }
        ]
        }
    ]

