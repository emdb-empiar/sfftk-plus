# -*- coding: utf-8 -*-
# parser.py
"""Parses command-line options"""
import argparse
import sys

from sfftk.core.parser import add_args
from sfftk.core.print_tools import print_date

__author__ = 'Paul K. Korir, PhD'
__email__ = 'pkorir@ebi.ac.uk, paul.korir@gmail.com'
__date__ = '2016-06-10'

verbosity_range = range(4)

Parser = argparse.ArgumentParser(prog='sffp', description="The EMDB-SFF Toolkit (sfftk)")

subparsers = Parser.add_subparsers(
    title='Tools',
    dest='subcommand',
    description='The EMDB-SFF Extended Toolkit (sfftk-plus) provides the following tools:',
    metavar="EMDB-SFF tools"
)

# ===============================================================================
# common arguments
# ===============================================================================
shipped_configs = {
    'args': ['-b', '--shipped-configs'],
    'kwargs': {
        'default': False,
        'action': 'store_true',
        'help': 'use shipped configs only if config path and user configs fail [default: False]'
    }
}
center = {
    'args': ['-c', '--center'],
    'kwargs': {
        'default': False,
        'action': 'store_true',
        'help': 'translate the segmentation to have its center around the origin [default: False]'
    }
}
dataset = {
    'args': ['-d', '--dataset'],
    'kwargs': {
        'default': None,
        'help': "dataset name [default: None]"
    }
}
details_param = {
    'args': ['-d', '--details'],
    'kwargs': {
        'default': "",
        'help': "populates <details>...</details> in the XML file [default: '']"
    }
}
# host = {
#     'args': ['-H', '--host'],
#     'kwargs': {
# #         'default': 'localhost',
#         'help': 'OMERO-server host'
#         }
#     }
image_id = {
    'args': ['-i', '--image-id'],
    'kwargs': {
        'type': int,
        'help': "id of image in OMERO-server"
    }
}
project = {
    'args': ['-j', '--project'],
    'kwargs': {
        'default': None,
        'help': "project name [default: None]"
    }
}
local = {
    'args': ['--local'],
    'kwargs': {
        'default': False,
        'action': 'store_true',
        'help': 'connect to local OMERO using credentials specified in configs [default: False; mutex with --remote]'
    }
}
mask_value = {
    'args': ['-m', '--mask-value'],
    'kwargs': {
        'default': 1,
        'type': int,
        'help': "for masks (threeDVolume segments): an integer value used mark masked regions [default: 1]"
    }
}
normals_off = {
    'args': ['-N', '--normals-off'],
    'kwargs': {
        'action': 'store_true',
        'default': False,
        'help': "do not to use normals if present [default: False]",
    }
}
image_name = {
    'args': ['-n', '--image-name'],
    'kwargs': {
        'default': None,
        'help': "name of the image of interest"
    }
}
output_path = {
    'args': ['-O', '--output-path'],
    'kwargs': {
        'default': './',
        'help': "path to which files will be written [default: ./]"
    }
}
output = {
    'args': ['-o', '--output'],
    'kwargs': {
        #         'type': argparse.FileType('w'),
        #         'default': sys.stdout,
        'required': True,
        'help': "file to convert to [default: sys.stdout]"
    }
}
# password = {
#     'args': ['-P', '--password'],
#     'kwargs': {
#         'help': 'OMERO-server password'
#         }
#     }
config_path = {
    'args': ['-p', '--config-path'],
    'kwargs': {
        'help': "path to configs file"
    }
}
# port = {
#     'args':['-Q', '--port'],
#     'kwargs': {
# #         'default': '4064',
#         'type': int,
#         'help': 'OMERO-server host'
#         }
#     }
primary_descriptor = {
    'args': ['-R', '--primary-descriptor'],
    'kwargs': {
        'default': None,
        'help': "populates the <primaryDescriptor>...</primaryDescriptor> to this value [valid values:  threeDVolume, contourList, meshList, shapePrimitiveList]"
    }
}
remote = {
    'args': ['--remote'],
    'kwargs': {
        'default': False,
        'action': 'store_true',
        'help': 'connect to remote OMERO using credentials specified in configs [default: False; mutex with --local]'
    }
}
summary = {
    'args': ['-s', '--summary'],
    'kwargs': {
        'default': False,
        'action': 'store_true',
        'help': "display a summary instead of details [default: False]"
    }
}
default_transparency = 0.5
transparency = {
    'args': ['-t', '--transparency'],
    'kwargs': {
        'type': float,
        'default': default_transparency,
        'help': "set the transparency on a scale of 0.0 (transparent) to 1.0 (opaque) [default: {}]".format(
            default_transparency)
    }
}
# user = {
#     'args': ['-U', '--user'],
#     'kwargs': {
#         'help': 'OMERO-server username',
#         }
#     }
verbose = {
    'args': ['-v', '--verbose'],
    'kwargs': {
        'action': 'store_true',
        'default': False,
        'help': "verbose output"
    }
}

# ===============================================================================
# updateschema subparser
# ===============================================================================
# updateschema_parser = subparsers.add_parser('updateschema', description="Update the schema API (schema/emdb_sff.py and schema/roi.py)", help="update schemas (emdb_sff.py or roi.py) using generateDS.py")
# updateschema_parser.add_argument('-s', '--schema-file', help="path to schema (.xsd) file")

# =============================================================================
# config subparser
# =============================================================================
config_parser = subparsers.add_parser(
    'config',
    description="Configuration utility",
    help="manage sfftkplus configs"
)

config_subparsers = config_parser.add_subparsers(
    title="sfftkplus configurations",
    dest="config_subcommand",
    description="Persistent configurations utility",
    metavar="Commands:"
)

# =============================================================================
# config: list
# =============================================================================
list_configs_parser = config_subparsers.add_parser(
    'list',
    description="List sfftkplus configuration parameters",
    help="list sfftkplus configs"
)
add_args(list_configs_parser, config_path)
add_args(list_configs_parser, shipped_configs)

# =============================================================================
# config: get
# =============================================================================
get_configs_parser = config_subparsers.add_parser(
    'get',
    description='Get the value of a single configuration parameter',
    help='get single sfftkplus config'
)
get_configs_parser.add_argument(
    'name', help="the name of the argument to retrieve"
)
add_args(get_configs_parser, config_path)
add_args(get_configs_parser, shipped_configs)

# =============================================================================
# config: set
# =============================================================================
set_configs_parser = config_subparsers.add_parser(
    'set',
    description='Set the value of a single configuration parameter',
    help='set single sfftkplus config'
)
set_configs_parser.add_argument(
    'name', help="the name of the argument to set",
)
set_configs_parser.add_argument(
    'value', help="the value of the argument to set",
)
add_args(set_configs_parser, config_path)
add_args(set_configs_parser, shipped_configs)

# =============================================================================
# config: del
# =============================================================================
del_configs_parser = config_subparsers.add_parser(
    'del',
    description='Delete the named configuration parameter',
    help='delete single sfftkplus config'
)
del_configs_parser.add_argument(
    'name', help="the name of the argument to be deleted"
)
add_args(del_configs_parser, config_path)
add_args(del_configs_parser, shipped_configs)

# =============================================================================
# config: clear
# =============================================================================
clear_configs_parser = config_subparsers.add_parser(
    'clear',
    description='Clear all configuration parameters',
    help='clear all sfftkplus configs'
)
add_args(clear_configs_parser, config_path)
add_args(clear_configs_parser, shipped_configs)

# ===============================================================================
# list subparser
# ===============================================================================
list_parser = subparsers.add_parser(
    'list',
    description="Lists various entities",
    help="list various entities"
)
list_type_group = list_parser.add_mutually_exclusive_group()
list_type_group.add_argument('-I', '--images', action='store_true', default=False, help="list images in OMERO-server")
list_type_group.add_argument('-R', '--rois', action='store_true', default=False, help="list ROIs in OMERO-server")
list_find_group = list_parser.add_mutually_exclusive_group()
add_args(list_find_group, image_id)
add_args(list_find_group, image_name)
add_args(list_parser, project)
add_args(list_parser, dataset)
add_args(list_parser, summary)
add_args(list_parser, config_path)
add_args(list_parser, shipped_configs)
add_args(list_parser, verbose)
# list_parser.add_argument(*user['args'], **user['kwargs'])
# list_parser.add_argument(*password['args'], **password['kwargs'])
# list_parser.add_argument(*host['args'], **host['kwargs'])
# list_parser.add_argument(*port['args'], **port['kwargs'])
local_or_remote = list_parser.add_mutually_exclusive_group()
add_args(local_or_remote, local)
add_args(local_or_remote, remote)
# group.add_argument(*local['args'], **local['kwargs'])
# group.add_argument(*remote['args'], **remote['kwargs'])


# ===============================================================================
# createroi subparser
# ===============================================================================
createroi_parser = subparsers.add_parser('createroi', description="Create ROIs and write to file",
                                         help="create ROIs and write to file")
createroi_parser.add_argument('sff_file', help="file containing segmentations to be converted into ROIs")
add_args(createroi_parser, config_path)
add_args(createroi_parser, shipped_configs)
createroi_parser.add_argument(*output['args'], **output['kwargs'])
createroi_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
createroi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
createroi_parser.add_argument(*transparency['args'], **transparency['kwargs'])
# createroi_parser.add_argument('-s', '--smooth', action='store_true', default=False, help="attempt to smoothen mesh [default: False]")
# createroi_parser.add_argument('-i', '--smooth-iterations', default=50, type=int, help="the number of smoothing iterations [default: 50]")
createroi_parser.add_argument('-I', '--image-name-root',
                              help="the root name of the file in OMERO; e.g. 'emd_1080-top.map' has image root 'emd_1080'")
createroi_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])
createroi_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])
createroi_parser.add_argument('-q', '--quick-pick', type=int,
                              help="if multiple IDs are found pick the one at the specified position; uses 1-based indexing for natural positioning [default: None]")

# ===============================================================================
# attachroi subparser
# ===============================================================================
attachroi_parser = subparsers.add_parser('attachroi',
                                         description="Associate the ROIs from the file with the image (by ID) in the OMERO-server",
                                         help="attach ROIs to image")
attachroi_parser.add_argument('roi_file', help="file with ROIs")
attachroi_parser.add_argument('-x', '--x-image-id', type=int, help="id of front image in OMERO-server")
attachroi_parser.add_argument('-y', '--y-image-id', type=int, help="id of rightside image in OMERO-server")
attachroi_parser.add_argument('-z', '--z-image-id', type=int, help="id of top image in OMERO-server")
# attachroi_parser.add_argument(*user['args'], **user['kwargs'])
# attachroi_parser.add_argument(*password['args'], **password['kwargs'])
# attachroi_parser.add_argument(*host['args'], **host['kwargs'])
# attachroi_parser.add_argument(*port['args'], **port['kwargs'])
attachroi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
attachroi_parser.add_argument(*config_path['args'], **config_path['kwargs'])
add_args(attachroi_parser, shipped_configs)
local_or_remote = attachroi_parser.add_mutually_exclusive_group()
add_args(local_or_remote, local)
add_args(local_or_remote, remote)

# ===============================================================================
# delroi subparser
# ===============================================================================
delroi_parser = subparsers.add_parser('delroi',
                                      description="Delete all ROIs associated with the image (by ID) in the OMERO-server",
                                      help="delete ROIs associated with image")
delroi_type_group = delroi_parser.add_mutually_exclusive_group()
delroi_type_group.add_argument('-r', '--roi-id', type=int, help="id of ROI in OMERO-server")
delroi_type_group.add_argument(*image_id['args'], **image_id['kwargs'])
# delroi_parser.add_argument(*user['args'], **user['kwargs'])
# delroi_parser.add_argument(*password['args'], **password['kwargs'])
# delroi_parser.add_argument(*host['args'], **host['kwargs'])
# delroi_parser.add_argument(*port['args'], **port['kwargs'])
delroi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
delroi_parser.add_argument(*config_path['args'], **config_path['kwargs'])
add_args(delroi_parser, shipped_configs)
local_or_remote = delroi_parser.add_mutually_exclusive_group()
add_args(local_or_remote, local)
add_args(local_or_remote, remote)

# ===============================================================================
# view3d subparser
# ===============================================================================
view3d_parser = subparsers.add_parser('view3d', description="View 3D rendering of segmentations",
                                      help="render 3D model")
view3d_parser.add_argument('sff_file', help="any SFF file")
view3d_parser.add_argument('-A', '--all-contours', action='store_true', default=False,
                           help="show all contours [default: False]")
view3d_parser.add_argument('-C', '--keep-contours', action='store_true', default=False,
                           help="do not convert contours to meshes; only applies for contourList [default: False]")
view3d_parser.add_argument('-F', '--full-screen', action='store_true', default=False,
                           help="show models in full-screen mode [default: False]")
view3d_parser.add_argument('-M', '--exclude-mesh', action='store_true', default=False,
                           help="do not display the main mesh [default: False]")
view3d_parser.add_argument('-X', '--x-contours', action='store_true', default=False,
                           help="show x contours [default: False]")
view3d_parser.add_argument('-Y', '--y-contours', action='store_true', default=False,
                           help="show y contours [default: False]")
view3d_parser.add_argument('-Z', '--z-contours', action='store_true', default=False,
                           help="show z contours [default: False]")
view3d_parser.add_argument('-a', '--no-orientation-axes', action='store_true', default=False,
                           help="do not display orientation axes (bottom right of viewport) [default: True]")
view3d_parser.add_argument('-B', '--background-colour', nargs=3, default=[0.1, 0.2, 0.4], type=float,
                           help="set the colour to be used for the background [default: 0.1 0.2, 0.4]")
view3d_parser.add_argument('-c', '--cube-axes', type=int,
                           help="how to display the cube axes; 0 - outer edges; 1 - closest triad; 2 - furthest triad; 3 - static triad; 4 - static edges; [default: None]")
view3d_parser.add_argument('-e', '--view-edges', action='store_true', default=False, help="show edges [default: False]")
view3d_parser.add_argument('-f', '--fill-holes', action='store_true', default=False,
                           help="attempt to fill holes in mesh [default: False]")
# view3d_parser.add_argument('-i', '--smooth-iterations', default=50, type=int, help="the number of smoothing iterations [default: 50]")
view3d_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])
# view3d_parser.add_argument('-s', '--smooth', action='store_true', default=False, help="attempt to smoothen mesh [default: False]")
view3d_parser.add_argument('-w', '--wireframe', action='store_true', default=False,
                           help="use wireframe representation (but retains contours as solids) [default: False]")
view3d_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
view3d_parser.add_argument(*transparency['args'], **transparency['kwargs'])
view3d_parser.add_argument(*verbose['args'], **verbose['kwargs'])
view3d_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])
add_args(view3d_parser, config_path)
add_args(view3d_parser, shipped_configs)

# ===============================================================================
# export
# ===============================================================================
export_parser = subparsers.add_parser('export', description="Export segmentation as various file formats",
                                      help="export as file")
export_parser.add_argument('sff_file', help="any SFF file")
export_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
export_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])
export_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])
export_parser.add_argument(*verbose['args'], **verbose['kwargs'])
export_parser.add_argument(*output_path['args'], **output_path['kwargs'])
add_args(export_parser, config_path)
add_args(export_parser, shipped_configs)
add_args(export_parser, transparency)
add_args(export_parser, center)

# get the full list of tools from the Parser object
tool_list = ['core', 'formats', 'omero', 'readers', 'schema', 'main']

# tests
test_help = "one or none of the following: {}".format(", ".join(tool_list))
tests_parser = subparsers.add_parser('tests', description="Run unit tests", help="run unit tests")
tests_parser.add_argument('tool', nargs='*', default='all', help=test_help)
tests_parser.add_argument('-v', '--verbosity', default=1, type=int,
                          help="set verbosity; valid values: %s [default: 0]" % ", ".join(map(str, verbosity_range)))
add_args(tests_parser, config_path)
add_args(tests_parser, shipped_configs)


# test_parser = subparsers.add_parser('test', description="Run unit tests", help="run unit tests")
# test_parser.add_argument('tool', nargs='*', default='all', help=test_help)
# test_parser.add_argument('-v', '--verbosity', default=1, type=int, help="set verbosity; valid values: %s [default: 0]" % ", ".join(map(str, verbosity_range)))


def parse_args(_args):
    """
    Parse and check command-line arguments
    
    Subcommand handlers defined in __main__.py (e.g. handle_conver(...)) should not have to check arguments for consistency
    
    :param list _args: list of arguments
    :return: parsed arguments
    :rtype: `argparse.Namespace`
    :return: config dict-like object
    :rtype: ``sfftk.core.configs.Config``
    """
    """
    :TODO: handle credentials in configs here instead of sffplus.py
    """
    # if we have no subcommands then show the available tools
    if len(_args) == 0:
        Parser.print_help()
        sys.exit(0)
    # if we only have a subcommand then show that subcommand's help
    elif len(_args) == 1:
        if _args[0] == "tests":
            pass
        elif _args[0] in Parser._actions[1].choices.keys():
            exec ('{}_parser.print_help()'.format(_args[0]))
            sys.exit(0)

    # parse args
    args = Parser.parse_args(_args)
    from sfftk.core.configs import load_configs
    from .configs import SFFPConfigs
    configs = load_configs(
        args,
        user_folder=".sfftkplus",
        conf_fn="sffp.conf",
        config_class=SFFPConfigs
    )
    if args.subcommand == 'list':
        # enforce local if specified
        if args.local:
            configs['CONNECT_WITH'] = 'LOCAL'
        elif args.remote:
            configs['CONNECT_WITH'] = 'REMOTE'

    # createroi
    elif args.subcommand == 'createroi':
        # ensure valid primary_descriptor
        if args.primary_descriptor:
            try:
                assert args.primary_descriptor in ['threeDVolume', 'contourList', 'meshList', 'shapePrimitiveList']
            except AssertionError:
                print_date('Invalid value for primaryDescriptor: %s' % args.primary_descriptor)
                return None, configs

        # ensure valid transparencey
        if args.transparency:
            try:
                assert 0 <= args.transparency <= 1
            except AssertionError:
                print_date("Invalid value for transparency: {}; should be between 0 and 1 (inclusive)".format(
                    args.transparency))
                return None, configs

        # ensure mask value is an integer (or long)
        if args.mask_value:
            try:
                assert isinstance(args.mask_value, int) or isinstance(args.mask_value, long)
            except AssertionError:
                print_date("Non-integer for mask value")
                return None, configs
        # quick pick values are 1-based (not 0-based)
        if args.quick_pick is not None:
            args.quick_pick -= 1
    # attachroi
    elif args.subcommand == 'attachroi':
        # enforce local if specified
        if args.local:
            configs['CONNECT_WITH'] = 'LOCAL'
        elif args.remote:
            configs['CONNECT_WITH'] = 'REMOTE'

    # delroi
    elif args.subcommand == 'delroi':
        # enforce local if specified
        if args.local:
            configs['CONNECT_WITH'] = 'LOCAL'
        elif args.remote:
            configs['CONNECT_WITH'] = 'REMOTE'

        # ensure that we have either an image or ROI ID
        if not args.image_id and not args.roi_id:
            raise ValueError('Missing either image (-i) or ROI (-r) ID')

        # ensure that both image and ROI ID are not set simultaneously
        if args.image_id and args.roi_id:
            raise ValueError('Only set one of image (-i) or ROI (-r) ID; not both')

    # tests
    elif args.subcommand == 'tests':
        if isinstance(args.tool, list):
            for tool in args.tool:
                try:
                    assert tool in tool_list
                except AssertionError:
                    print >> sys.stderr, "Unknown tool: {}".format(tool)
                    print >> sys.stderr, "Available tools for test: {}".format(", ".join(tool_list))

        if args.verbosity:
            try:
                assert args.verbosity in range(4)
            except:
                raise ValueError("Verbosity should be in %s-%s: %s given" % (
                verbosity_range[0], verbosity_range[-1], args.verbosity))

    # view3d
    elif args.subcommand == 'view3d':
        if args.primary_descriptor:
            try:
                assert args.primary_descriptor in ['threeDVolume', 'contourList', 'meshList', 'shapePrimitiveList']
            except:
                raise ValueError('Invalid value for primaryDescriptor: %s' % args.primary_descriptor)

        # ensure valid transparencey
        assert 0 <= args.transparency <= 1

        #  ensure valid background colours
        assert 0 <= args.background_colour[0] <= 1
        assert 0 <= args.background_colour[1] <= 1
        assert 0 <= args.background_colour[2] <= 1

        # view contours
        # don't specify -A and -X, -Y, and/or -Z
        if args.all_contours or args.x_contours or args.y_contours or args.z_contours:
            assert (args.all_contours and not (args.x_contours or args.y_contours or args.z_contours)) or (
            not args.all_contours and (args.x_contours or args.y_contours or args.z_contours))

        # cube axes validity
        assert (0 <= args.cube_axes <= 4) or (args.cube_axes is None)

    return args, configs
