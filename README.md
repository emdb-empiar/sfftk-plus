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