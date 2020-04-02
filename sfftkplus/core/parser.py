# -*- coding: utf-8 -*-
# parser.py
"""Parses command-line options"""
from __future__ import print_function

import os

from sfftk.core import parser # convert_parser, prep_parser, config_parser, notes_parser
from sfftkrw.core import _dict_iter_keys
from sfftkrw.core.parser import add_args
from sfftkrw.core.print_tools import print_date

__author__ = 'Paul K. Korir, PhD'
__email__ = 'pkorir@ebi.ac.uk, paul.korir@gmail.com'
__date__ = '2016-06-10'

verbosity_range = range(4)

parser.Parser.description = u"The Extended EMDB-SFF Toolkit (sfftk-plus)"

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
# details_param = {
#     'args': ['-d', '--details'],
#     'kwargs': {
#         'default': "",
#         'help': "populates <details>...</details> in the XML file [default: '']"
#     }
# }
FORMAT_LIST = [
    ('roi', 'ROI'),
    ('json', 'JSON'),
]
format_ = {
    'args': ['-f', '--format'],
    'kwargs': {
        'default': FORMAT_LIST[0][0],
        'choices': list(map(lambda f: f[0], FORMAT_LIST)),
        'help': "output file format; valid options are: {} [default: roi]".format(
            ", ".join(map(lambda x: "{} ({})".format(x[0], x[1]), FORMAT_LIST))
        ),
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
        'help': "for masks (three_d_volume segments): an integer value used mark masked regions [default: 1]"
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
        # 'required': True,
        'help': "the name of the output file; by default this is the name of the input file with an .roi extension",
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
        'help': "populates the <primary_descriptor>...</primary_descriptor> to this value [valid values:  three_d_volume, mesh_list, shape_primitive_list]"
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
default_transparency = 1.0
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
overwrite = {
    'args': ['-w', '--overwrite'],
    'kwargs': {
        'default': False,
        'action': 'store_true',
        'help': "overwrite the output [default: False]",
    }
}

# ===============================================================================
# list subparser
# ===============================================================================
list_parser = parser.subparsers.add_parser(
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
roi_parser = parser.subparsers.add_parser('roi',
                                   description='Work with 2D regions-of-interest (ROIs). You can create, attach and delete ROIs to an OMERO server',
                                   help='work with 2D ROIs')

roi_subparser = roi_parser.add_subparsers(
    title='Region-of-Interest (ROI) tools',
    dest='roi_subcommand',
    description='Convert 3D segmentations to 2D ROIs associated with each slice of the corresponding image data',
    metavar='ROI tools',
)

create_roi_parser = roi_subparser.add_parser('create', description="Create ROIs and write to file",
                                             help="create ROIs and write to file")
add_args(create_roi_parser, config_path)
add_args(create_roi_parser, shipped_configs)
create_roi_parser.add_argument(*output['args'], **output['kwargs'])
add_args(create_roi_parser, overwrite)
add_args(create_roi_parser, format_)
create_roi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
# mutex parser group
image_name_root_or_xyz_createroi_parser = create_roi_parser.add_mutually_exclusive_group()
image_name_root_or_xyz_createroi_parser.add_argument('-I', '--image-name-root',
                                                     help="the root name of the file in OMERO; e.g. 'emd_1080-top.map' has image root 'emd_1080'")
image_name_root_or_xyz_createroi_parser.add_argument(
    '--top-front-right',
    nargs=3,
    type=int,
    metavar=('TOP-IMAGE-ID', 'FRONT-IMAGE-ID', 'RIGHT-IMAGE-ID'),
    help="explicit image IDs for top, front and right (side) perspectives, respectively"
)
create_roi_parser.add_argument(
    '-q', '--quick-pick',
    type=int,
    default=1,
    help="if multiple IDs are found pick the one at the specified position; "
         "uses 1-based indexing for natural positioning [default: 1]"
)
create_roi_parser.add_argument(
    'sff_file',
    help="file containing segmentations to be converted into ROIs; this could also be an ROI file (*.roi) for "
         "modifying the image ids"
)
create_roi_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
create_roi_parser.add_argument(*transparency['args'], **transparency['kwargs'])
create_roi_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])
create_roi_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])

create_roi_parser.add_argument(
    '-i', '--reset-ids',
    action='store_true',
    default=False,
    help="modify the image IDs based on a fresh search [default: False]"
)

# ===============================================================================
# attachroi subparser
# ===============================================================================
attach_roi_parser = roi_subparser.add_parser('attach',
                                             description="Associate the ROIs from the file with the image (by ID) in the OMERO-server",
                                             help="attach ROIs to image")
attach_roi_parser.add_argument('roi_file', help="file with ROIs")
attach_roi_parser.add_argument('-x', '--x-image-id', type=int, help="id of front image in OMERO-server")
attach_roi_parser.add_argument('-y', '--y-image-id', type=int, help="id of rightside image in OMERO-server")
attach_roi_parser.add_argument('-z', '--z-image-id', type=int, help="id of top image in OMERO-server")
# attachroi_parser.add_argument(*user['args'], **user['kwargs'])
# attachroi_parser.add_argument(*password['args'], **password['kwargs'])
# attachroi_parser.add_argument(*host['args'], **host['kwargs'])
# attachroi_parser.add_argument(*port['args'], **port['kwargs'])
attach_roi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
attach_roi_parser.add_argument(*config_path['args'], **config_path['kwargs'])
add_args(attach_roi_parser, shipped_configs)
local_or_remote = attach_roi_parser.add_mutually_exclusive_group()
add_args(local_or_remote, local)
add_args(local_or_remote, remote)

# ===============================================================================
# delroi subparser
# ===============================================================================
del_roi_parser = roi_subparser.add_parser('del',
                                          description="Delete all ROIs associated with the image (by ID) in the OMERO-server",
                                          help="delete ROIs associated with image")
delroi_type_group = del_roi_parser.add_mutually_exclusive_group()
delroi_type_group.add_argument('-r', '--roi-id', type=int, help="id of ROI in OMERO-server")
delroi_type_group.add_argument(*image_id['args'], **image_id['kwargs'])
# delroi_parser.add_argument(*user['args'], **user['kwargs'])
# delroi_parser.add_argument(*password['args'], **password['kwargs'])
# delroi_parser.add_argument(*host['args'], **host['kwargs'])
# delroi_parser.add_argument(*port['args'], **port['kwargs'])
del_roi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
del_roi_parser.add_argument(*config_path['args'], **config_path['kwargs'])
add_args(del_roi_parser, shipped_configs)
local_or_remote = del_roi_parser.add_mutually_exclusive_group()
add_args(local_or_remote, local)
add_args(local_or_remote, remote)

# ===============================================================================
# view3d subparser
# ===============================================================================
parser.view_parser.add_argument('--visualise', action='store_true', default=False, help='display a 3D rendering of the '
                                                                                 'segmentation geometry in the provided '
                                                                                 'file')
parser.view_parser.add_argument('-A', '--all-contours', action='store_true', default=False,
                         help="show all contours [default: False]")
parser.view_parser.add_argument('-F', '--full-screen', action='store_true', default=False,
                         help="show models in full-screen mode [default: False]")
parser.view_parser.add_argument('-M', '--exclude-mesh', action='store_true', default=False,
                         help="do not display the main mesh [default: False]")
parser.view_parser.add_argument('-X', '--x-contours', action='store_true', default=False,
                         help="show x contours [default: False]")
parser.view_parser.add_argument('-Y', '--y-contours', action='store_true', default=False,
                         help="show y contours [default: False]")
parser.view_parser.add_argument('-Z', '--z-contours', action='store_true', default=False,
                         help="show z contours [default: False]")
parser.view_parser.add_argument('-a', '--no-orientation-axes', action='store_true', default=False,
                         help="do not display orientation axes (bottom right of viewport) [default: True]")
parser.view_parser.add_argument('-B', '--background-colour', nargs=3, default=[0.1, 0.2, 0.4], type=float,
                         help="set the colour to be used for the background [default: 0.1 0.2, 0.4]")
parser.view_parser.add_argument('-c', '--cube-axes', type=int,
                         help="how to display the cube axes; 0 - outer edges; 1 - closest triad; 2 - furthest triad; 3 - static triad; 4 - static edges; [default: None]")
parser.view_parser.add_argument('-e', '--view-edges', action='store_true', default=False, help="show edges [default: False]")
parser.view_parser.add_argument('-f', '--fill-holes', action='store_true', default=False,
                         help="attempt to fill holes in mesh [default: False]")
parser.view_parser.add_argument('-i', '--smooth-iterations', default=50, type=int,
                         help="the number of smoothing iterations [default: 50]")
parser.view_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])
parser.view_parser.add_argument('-s', '--smooth', action='store_true', default=False,
                         help="attempt to smoothen mesh [default: False]")
parser.view_parser.add_argument('-w', '--wireframe', action='store_true', default=False,
                         help="use wireframe representation (but retains contours as solids) [default: False]")
parser.view_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
parser.view_parser.add_argument(*transparency['args'], **transparency['kwargs'])
parser.view_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])

# ===============================================================================
# export
# ===============================================================================
export_parser = parser.subparsers.add_parser('export', description="Export segmentation as various file formats",
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
parser.tool_list += ['all_sfftk_plus', 'formats_sfftk_plus', 'omero', 'readers_sfftk_plus', 'schema_sfftk_plus', 'main_sfftk_plus']

# tests
test_help = "one or none of the following: {}".format(", ".join(parser.tool_list))
tests_parser_tools = parser.tests_parser._actions[1]
tests_parser_tools.help = test_help


def parse_args(_args, use_shlex=False):
    """
    Parse and check command-line arguments
    
    Subcommand handlers defined in __main__.py (e.g. handle_conver(...)) should not have to check arguments for consistency
    
    :param str _args: command string
    :return: parsed arguments
    :rtype: `argparse.Namespace`
    :return: config dict-like object
    :rtype: ``sfftk.core.configs.Configs``
    """
    # use shlex
    if use_shlex:
        try:
            assert isinstance(_args, str)
        except AssertionError:
            return os.EX_USAGE, None
        import shlex
        _args = shlex.split(_args)

    if len(_args) > 0:
        if _args[0] in ['convert', 'notes', 'prep', 'config']:
            return parser.parse_args(_args)

    """
    :TODO: handle credentials in configs here instead of sffplus.py
    """
    # if we have no subcommands then show the available tools
    if len(_args) == 0:
        parser.Parser.print_help()
        return os.EX_OK, None
    # if we only have a subcommand then show that subcommand's help
    elif len(_args) == 1:
        if _args[0] == '-V' or _args[0] == '--version':
            from .. import SFFTKPLUS_VERSION
            print_date("sfftk-plus version: {}".format(SFFTKPLUS_VERSION))
            return os.EX_OK, None
        elif _args[0] in parser.Parser._actions[2].choices.keys():
            try:
                exec('parser.{}_parser.print_help()'.format(_args[0]))
            except AttributeError:
                exec('{}_parser.print_help()'.format(_args[0]))
            return os.EX_OK, None
    elif len(_args) == 2:
        if _args[0] == 'roi':
            if _args[1] == 'create':
                if _args[1] in _dict_iter_keys(parser.Parser._actions[2].choices['roi']._actions[1].choices):
                    exec('{}_roi_parser.print_help()'.format(_args[1]))
                    return os.EX_OK, None
            elif _args[1] == 'attach':
                if _args[1] in _dict_iter_keys(parser.Parser._actions[2].choices['roi']._actions[1].choices):
                    exec('{}_roi_parser.print_help()'.format(_args[1]))
                    return os.EX_OK, None
            elif _args[0] == 'del':
                if _args[1] in _dict_iter_keys(parser.Parser._actions[2].choices['roi']._actions[1].choices):
                    exec('{}_roi_parser.print_help()'.format(_args[1]))
                    return os.EX_OK, None
        elif _args[0] == 'notes':
            if _args[1] in _dict_iter_keys(parser.Parser._actions[2].choices['notes']._actions[1].choices):
                exec('parser.{}_notes_parser.print_help()'.format(_args[1]))
                return os.EX_OK, None
        elif _args[0] == 'prep':
            if _args[1] in _dict_iter_keys(parser.Parser._actions[2].choices['prep']._actions[1].choices):
                exec('parser.{}_prep_parser.print_help()'.format(_args[1]))
                return os.EX_OK, None
        elif _args[0] == 'config':
            if _args[1] in _dict_iter_keys(parser.Parser._actions[2].choices['config']._actions[1].choices):
                exec('parser.{}_config_parser.print_help()'.format(_args[1]))
                return os.EX_OK, None

    # parse args
    args = parser.Parser.parse_args(_args)
    from sfftk.core.configs import get_config_file_path, load_configs, Configs
    config_file_path = get_config_file_path(args, user_conf_fn='sff.conf', user_folder='~/.sfftk', config_class=Configs)
    configs = load_configs(config_file_path, config_class=Configs)
    # config
    if args.subcommand == 'config':
        if args.verbose:
            print_date("Reading configs from {}...".format(config_file_path))
            # handle config-specific argument modifications here
        if args.config_subcommand == 'del':
            if args.name not in configs:
                print_date("Missing config with name '{}'. Aborting...".format(args.name))
                return os.EX_USAGE, configs
            # if force pass
            if not args.force:
                default_choice = 'n'
                # get user choice
                user_choice = input("Are you sure you want to delete config '{}' [y/N]? ".format(
                    args.name)).lower()
                if user_choice == '':
                    choice = default_choice
                elif user_choice == 'n' or user_choice == 'N':
                    choice = 'n'
                elif user_choice == 'y' or user_choice == 'Y':
                    choice = 'y'
                else:
                    print_date("Invalid choice: '{}'")
                    return os.EX_DATAERR, configs
                # act on user choice
                if choice == 'n':
                    print_date("You have opted to cancel deletion of '{}'".format(args.name))
                    return os.EX_OK, configs
                elif choice == 'y':
                    pass
        elif args.config_subcommand == 'set':
            if args.name in configs:
                # if force pass
                if not args.force:
                    default_choice = 'n'
                    # get user choice
                    user_choice = input("Are you sure you want to overwrite config '{}={}' [y/N]? ".format(
                        args.name, configs[args.name])).lower()
                    if user_choice == '':
                        choice = default_choice
                    elif user_choice == 'n' or user_choice == 'N':
                        choice = 'n'
                    elif user_choice == 'y' or user_choice == 'Y':
                        choice = 'y'
                    else:
                        print_date("Invalid choice: '{}'")
                        return os.EX_DATAERR, configs
                    # act on user choice
                    if choice == 'n':
                        print_date("You have opted to cancel overwriting of '{}'".format(args.name))
                        return os.EX_OK, configs
                    elif choice == 'y':
                        pass
    elif args.subcommand == 'list':
        # enforce local if specified
        if args.local:
            configs['CONNECT_WITH'] = 'LOCAL'
        elif args.remote:
            configs['CONNECT_WITH'] = 'REMOTE'

    # roi
    elif args.subcommand == 'roi':
        # create
        if args.roi_subcommand == 'create':
            # make the output file name and check if it exists
            if args.output is None:
                ofn = '.'.join(args.sff_file.split('.')[:-1]) + '.{}'.format(args.format)
                if os.path.exists(ofn) and not args.overwrite:
                    print_date("Output file exists. Use --overwrite to replace it.")
                    return os.EX_USAGE, configs
                else:
                    print_date("Using output file {}".format(ofn))
                    args.output = ofn
            # ensure valid primary_descriptor
            if args.primary_descriptor:
                try:
                    assert args.primary_descriptor in ['three_d_volume', 'mesh_list', 'shape_primitive_list']
                except AssertionError:
                    print_date('Invalid value for primary_descriptor: %s' % args.primary_descriptor)
                    return os.EX_DATAERR, configs

            # ensure valid transparencey
            if args.transparency:
                try:
                    assert 0 <= args.transparency <= 1
                except AssertionError:
                    print_date("Invalid value for transparency: {}; should be between 0 and 1 (inclusive)".format(
                        args.transparency))
                    return os.EX_DATAERR, configs

            # ensure mask value is an integer (or long)
            if args.mask_value:
                try:
                    assert isinstance(args.mask_value, int)  # or isinstance(args.mask_value, long)
                except AssertionError:
                    print_date("Non-integer for mask value")
                    return os.EX_DATAERR, configs

            # quick pick values are 1-based (not 0-based)
            # if args.quick_pick is not None:
            if args.quick_pick <= 0:
                print_date("Invalid value for --quick-pick. Should be 1-based value of item in list e.g. the value of "
                           "'a' in ['a', 'b'] is 1 (one).")
                return os.EX_DATAERR, configs
            else:
                args.quick_pick -= 1  # make it a 0-based index for internal use

            # if we don't have --top-front-right set
            # then we can choose a default for -I/--image-name-root
            if not args.top_front_right:
                if args.image_name_root is None:
                    image_name_root = os.path.basename('.'.join(args.sff_file.split('.')[:-1]))
                    args.image_name_root = image_name_root
                    print_date("Setting image name root to {}".format(image_name_root))

        # attach
        elif args.subcommand == 'attach':
            # enforce local if specified
            if args.local:
                configs['CONNECT_WITH'] = 'LOCAL'
            elif args.remote:
                configs['CONNECT_WITH'] = 'REMOTE'

        # del
        elif args.subcommand == 'del':
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
        # normalise tool list
        # if 'all' is specified together with others then it should simply be 'all'
        if 'all' in args.tool:
            args.tool = ['all']
        for tool in args.tool:
            try:
                assert tool in parser.tool_list
            except AssertionError:
                print_date(
                    "Unknown tool: {}; Available tools for test: {}".format(tool, ", ".join(parser.tool_list))
                )
                return os.EX_USAGE, configs
        if args.verbosity:
            try:
                assert args.verbosity in range(4)
            except:
                print_date(
                    "Verbosity should be in {}-{}: {} given".format(
                        verbosity_range[0],
                        verbosity_range[-1],
                        args.verbosity
                    )
                )
                return os.EX_USAGE, configs
        # if args.verbosity:
        #     try:
        #         assert args.verbosity in range(4)
        #     except:
        #         raise ValueError("Verbosity should be in %s-%s: %s given" % (
        #             verbosity_range[0], verbosity_range[-1], args.verbosity))

    # view3d
    elif args.subcommand == 'view':
        if args.primary_descriptor:
            try:
                assert args.primary_descriptor in ['three_d_volume', 'mesh_list', 'shape_primitive_list']
            except:
                raise ValueError('Invalid value for primary_descriptor: %s' % args.primary_descriptor)

        # ensure valid transparency
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
        if args.cube_axes is not None:
            assert (0 <= args.cube_axes <= 4)

    return args, configs
