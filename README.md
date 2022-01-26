# `sfftk-plus`
Extended toolkit for working with EMDB-SFF files

## Outline

- What does this package actually do?
    - Create ROI as:
        - XML files
        - JSON files
    - Export segmentations as VTP
    - Attach ROI to images in OMERO instance
- How should we set it up to work correctly?
    - Install `python-omero`
    - Install `python-vtk`

## Installation

You can install some of the dependencies using Anaconda
These will depend on the OMERO version because clients and server versions must match.

For example

```bash
conda install -c bioconda zeroc-ice==3.6 python-omero==5.3.3  
```

## Add Configs

Your installation of `sfftk-plus` will need to have some settings in place before it can work correctly.

```bash
sffp config get --all
```

will list the available configs.

```bash
CONNECT_WITH         = REMOTE
OMERO_LOCAL_HOST     = localhost
OMERO_LOCAL_USER     =
OMERO_LOCAL_PASSWORD =
OMERO_LOCAL_PORT     = 4064
OMERO_REMOTE_HOST    = <remote_omero_host>
OMERO_REMOTE_USER    = <remote_omero_user>
OMERO_REMOTE_PASSWORD = <remote_omero_password>
OMERO_REMOTE_PORT    = 4064
IMAGE_DB_LOCAL_NAME  =
IMAGE_DB_LOCAL_HOST  =
IMAGE_DB_LOCAL_USER  =
IMAGE_DB_LOCAL_PASS  =
IMAGE_DB_LOCAL_PORT  = 5432
IMAGE_DB_REMOTE_NAME = <remote_pg_db_name>
IMAGE_DB_REMOTE_HOST = <remote_pg_host>
IMAGE_DB_REMOTE_USER = <remote_pg_user>
IMAGE_DB_REMOTE_PASS = <remote_pg_password>
IMAGE_DB_REMOTE_PORT = 5432
```

Please set the configs as follows:

```bash
sff config set OMERO_REMOTE_HOST root
```

## Using `sfftk-plus`

> These notes are part of an email that I sent to someone who was interested in using `sfftk-plus`. They are likely to be incomplete and will be updated as and when the need arises.

First, I recommend installing sfftk-plus into a virtual environment of your choice.

Install it with `pip install git+https://github.com/emdb-empiar/sfftk-plus`.

This will also install the base libraries `sfftk-rw` and `sfftk` in addition to dependencies `vtk` and others. `vtk` is the core graphics library used to do geometrical manipulation.

`sfftk-plus` uses commands on the CLI. The main entrypoint is `sff`. This should display something like this:

```shell
usage: sff [-h] [-V] EMDB-SFF Read/Write Tools ...

The Extended EMDB-SFF Toolkit (sfftk-plus)

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show the sfftk-rw version string and the supported
                        EMDB-SFF Read/Write version string

Tools:
  The EMDB-SFF Read/Write Toolkit (sff) provides the following tools:

  EMDB-SFF Read/Write Tools
    convert             converts between EMDB-SFF formats
    view                view file summary
    tests               run unit tests
    prep                prepares a segmentation
    config              manage sfftk configs
    notes               annotate an EMDB-SFF file
    list                list various entities
    roi                 work with 2D ROIs
    export              export as file
```

If you are interested in this package then it is likely that you need to use the `view`, `roi` and `export` commands. Other commands may be used as described in `sfftk-rw` (https://sfftk-rw.readthedocs.io/en/latest/index.html) and `sfftk` (https://sfftk.readthedocs.io/en/latest/) manuals on EMDB-SFF files. Also keep in mind that both `sfftk-rw` and `sfftk` provide a programmatic API to create EMDB-SFF files from scratch!

You can use the `view` command to render the segmentation in 3D (use the `--visualise` option). Run `sff view` to display the complete help options.

You can use the `sff roi create` to create ROI slices. It requires an EMDB-SFF file. Other options are optional that we use in our backend. Use the `-f json -o emd_1547.json` option to get ROIs in JSON. The data in those files can be rendered directly using `raphael.js` (https://dmitrybaranovskiy.github.io/raphael/) with the appropriate SVG container markup.

You can use `sff export` to export VTP files. VTP files are native VTK files that capture compact geometry. Read more about VTP files at https://kitware.github.io/vtk-js/examples/GeometryViewer.html.
