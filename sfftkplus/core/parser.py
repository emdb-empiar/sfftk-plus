# -*- coding: utf-8 -*-
# parser.py
"""Parses command-line options"""
import argparse
import sys

from sfftk.core.print_tools import print_date


__author__  = 'Paul K. Korir, PhD'
__email__   = 'pkorir@ebi.ac.uk, paul.korir@gmail.com'
__date__    = '2016-06-10'




verbosity_range = range(4)

Parser = argparse.ArgumentParser(prog='sffp', description="The EMDB-SFF Toolkit (sfftk)")

subparsers = Parser.add_subparsers(
    title='Tools', 
    dest='subcommand', 
    description='The EMDB-SFF Extended Toolkit (sfftk-plus) provides the following tools:', 
    metavar="EMDB-SFF tools"
    )

#===============================================================================
# common arguments
#===============================================================================
use_local = {
    'args': ['-l', '--local'],
    'kwargs': {
        'default': True,
        'action': 'store_true',
        'help': 'connect to local OMERO using credentials specified in configs [default: True; mutex with -r/--remote]'
        }
    }
use_remote = {
    'args': ['-r', '--remote'],
    'kwargs': {
        'default': False,
        'action': 'store_true',
        'help': 'connect to remote OMERO using credentials specified in configs [default: False; mutex with -l/--local]'
        }
    }
user = {
    'args': ['-U', '--user'], 
    'kwargs': {
        'help': 'OMERO-server username',
        }
    }
password = {
    'args': ['-P', '--password'],
    'kwargs': {
        'help': 'OMERO-server password'
        }
    }
host = {
    'args': ['-H', '--host'],
    'kwargs': {
#         'default': 'localhost', 
        'help': 'OMERO-server host'
        }
    }
port = {
    'args':['-p', '--port'],
    'kwargs': {
#         'default': '4064', 
        'type': int,
        'help': 'OMERO-server host'
        }
    }
image_id = {
    'args':['-i', '--image-id'],
    'kwargs': {
        'type': int, 
        'help': "id of image in OMERO-server"
        }
    }
image_name = {
    'args': ['-n', '--image-name'],
    'kwargs': {
        'default': None,
        'help': "name of the image of interest"
        }
    }
primary_descriptor   = {
    'args':['-P', '--primary-descriptor'],
    'kwargs': {
        'default': None,
        'help': "populates the <primaryDescriptor>...</primaryDescriptor> to this value [valid values:  threeDVolume, contourList, meshList, shapePrimitiveList]"
        }
    }
output = {
    'args':['-o', '--output'],
    'kwargs': {
#         'type': argparse.FileType('w'),
#         'default': sys.stdout,
        'required': True,
        'help': "file to convert to [default: sys.stdout]"
        }
    }

output_path = {
    'args': ['-O', '--output-path'],
    'kwargs': {
        'default': './',
        'help': "path to which files will be written [default: ./]"
        }
    }

verbose = {
    'args':['-v', '--verbose'],
    'kwargs': {
        'action': 'store_true', 
        'default': False, 
        'help': "verbose output"
        }
    }
config_path = {
    'args':['-c', '--config-path'],
    'kwargs': {
        'help': "path to configs file"
        }
    }
details_param = {
    'args':['-d', '--details'],
    'kwargs': {
        'default': "", 
        'help': "populates <details>...</details> in the XML file [default: '']"
        }
    }
transparency = {
    'args':['-t', '--transparency'],
    'kwargs': {
        'type': float, 
        'default': 0.5, 
        'help': "set the transparency on a scale of 0.0 (transparent) to 1.0 (opaque) [default: 1.0]"
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
    'args': ['-n', '--normals-off'],
    'kwargs': {
        'action': 'store_true', 
        'default': False, 
        'help': "do not to use normals if present [default: False]",
        } 
    }

#===============================================================================
# updateschema subparser
#===============================================================================
updateschema_parser = subparsers.add_parser('updateschema', description="Update the schema API (schema/emdb_sff.py and schema/roi.py)", help="update schemas (emdb_sff.py or roi.py) using generateDS.py")
updateschema_parser.add_argument('-s', '--schema-file', help="path to schema (.xsd) file")

#===============================================================================
# list subparser
#===============================================================================
list_parser = subparsers.add_parser('list', description="Lists various entities", help="list various entities")
list_parser.add_argument('-I', '--images', action='store_true', default=False, help="list images in OMERO-server")
list_parser.add_argument('-R', '--rois', action='store_true', default=False, help="list ROIs in OMERO-server")
list_parser.add_argument(*image_id['args'], **image_id['kwargs'])
list_parser.add_argument(*image_name['args'], **image_name['kwargs'])
list_parser.add_argument(*user['args'], **user['kwargs'])
list_parser.add_argument(*password['args'], **password['kwargs'])
list_parser.add_argument(*host['args'], **host['kwargs'])
list_parser.add_argument(*port['args'], **port['kwargs'])
list_parser.add_argument(*verbose['args'], **verbose['kwargs'])
list_parser.add_argument(*config_path['args'], **config_path['kwargs'])
list_parser.add_argument('-j', '--project', default=None, help="project name [default: None]")
list_parser.add_argument('-d', '--dataset', default=None, help="dataset name [default: None]")
list_parser.add_argument('-s', '--summary', default=False, action='store_true', help="display a summary; not details [default: False]")
group = list_parser.add_mutually_exclusive_group()
group.add_argument(*use_local['args'], **use_local['kwargs'])
group.add_argument(*use_remote['args'], **use_remote['kwargs'])
# list supported file formats
# list configs

#===============================================================================
# createroi subparser
#===============================================================================
createroi_parser = subparsers.add_parser('createroi', description="Create ROIs and write to file", help="create ROIs and write to file")
createroi_parser.add_argument('sff_file', help="file containing segmentations to be converted into ROIs")
createroi_parser.add_argument(*output['args'], **output['kwargs'])
createroi_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
createroi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
createroi_parser.add_argument(*transparency['args'], **transparency['kwargs'])
# createroi_parser.add_argument('-s', '--smooth', action='store_true', default=False, help="attempt to smoothen mesh [default: False]")
# createroi_parser.add_argument('-i', '--smooth-iterations', default=50, type=int, help="the number of smoothing iterations [default: 50]")
createroi_parser.add_argument('-I', '--image-name-root', help="the root name of the file in OMERO; e.g. 'emd_1080-top.map' has image root 'emd_1080'")
createroi_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])
createroi_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])

#===============================================================================
# attachroi subparser
#===============================================================================
attachroi_parser = subparsers.add_parser('attachroi',  description="Associate the ROIs from the file with the image (by ID) in the OMERO-server", help="attach ROIs to image")
attachroi_parser.add_argument('roi_file', help="file with ROIs")
attachroi_parser.add_argument('-x', '--x-image-id', type=int, help="id of front image in OMERO-server")
attachroi_parser.add_argument('-y', '--y-image-id', type=int, help="id of rightside image in OMERO-server")
attachroi_parser.add_argument('-z', '--z-image-id', type=int, help="id of top image in OMERO-server")
attachroi_parser.add_argument(*user['args'], **user['kwargs'])
attachroi_parser.add_argument(*password['args'], **password['kwargs'])
attachroi_parser.add_argument(*host['args'], **host['kwargs'])
attachroi_parser.add_argument(*port['args'], **port['kwargs'])
attachroi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
attachroi_parser.add_argument(*config_path['args'], **config_path['kwargs'])

#===============================================================================
# delroi subparser
#===============================================================================
delroi_parser = subparsers.add_parser('delroi', description="Delete all ROIs associated with the image (by ID) in the OMERO-server", help="delete ROIs associated with image")
delroi_parser.add_argument('-r', '--roi-id', type=int, help="id of ROI in OMERO-server")
delroi_parser.add_argument(*image_id['args'], **image_id['kwargs'])
delroi_parser.add_argument(*user['args'], **user['kwargs'])
delroi_parser.add_argument(*password['args'], **password['kwargs'])
delroi_parser.add_argument(*host['args'], **host['kwargs'])
delroi_parser.add_argument(*port['args'], **port['kwargs'])
delroi_parser.add_argument(*verbose['args'], **verbose['kwargs'])
delroi_parser.add_argument(*config_path['args'], **config_path['kwargs'])

#===============================================================================
# config subparser
#===============================================================================
# config_parser = subparsers.add_parser('config', description="Write configs used by updateschema and OMERO-server connection", help="manage configs")
# config_parser.add_argument('config_name', help="name of config to write")
# config_parser.add_argument('config_value', help="value of config to write")

#===============================================================================
# view3d subparser
#===============================================================================
view3d_parser = subparsers.add_parser('view3d', description="View 3D rendering of segmentations", help="render 3D model")
view3d_parser.add_argument('sff_file', help="any SFF file")
view3d_parser.add_argument('-A', '--all-contours', action='store_true', default=False, help="show all contours [default: False]")
view3d_parser.add_argument('-C', '--keep-contours', action='store_true', default=False, help="do not convert contours to meshes; only applies for contourList [default: False]")
view3d_parser.add_argument('-F', '--full-screen', action='store_true', default=False, help="show models in full-screen mode [default: False]")
view3d_parser.add_argument('-M', '--exclude-mesh', action='store_true', default=False, help="do not display the main mesh [default: False]")
view3d_parser.add_argument('-X', '--x-contours', action='store_true', default=False, help="show x contours [default: False]")
view3d_parser.add_argument('-Y', '--y-contours', action='store_true', default=False, help="show y contours [default: False]")
view3d_parser.add_argument('-Z', '--z-contours', action='store_true', default=False, help="show z contours [default: False]")
view3d_parser.add_argument('-a', '--no-orientation-axes', action='store_true', default=False, help="do not display orientation axes (bottom right of viewport) [default: True]")
view3d_parser.add_argument('-b', '--background-colour', nargs=3, default=[0.1, 0.2, 0.4], type=float, help="set the colour to be used for the background [default: 0.1 0.2, 0.4]")
view3d_parser.add_argument('-c', '--cube-axes', type=int, help="how to display the cube axes; 0 - outer edges; 1 - closest triad; 2 - furthest triad; 3 - static triad; 4 - static edges; [default: None]")
view3d_parser.add_argument('-e', '--view-edges', action='store_true', default=False, help="show edges [default: False]")
view3d_parser.add_argument('-f', '--fill-holes', action='store_true', default=False, help="attempt to fill holes in mesh [default: False]")
# view3d_parser.add_argument('-i', '--smooth-iterations', default=50, type=int, help="the number of smoothing iterations [default: 50]")
view3d_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])
# view3d_parser.add_argument('-s', '--smooth', action='store_true', default=False, help="attempt to smoothen mesh [default: False]")
view3d_parser.add_argument('-w', '--wireframe', action='store_true', default=False, help="use wireframe representation (but retains contours as solids) [default: False]")
view3d_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
view3d_parser.add_argument(*transparency['args'], **transparency['kwargs'])
view3d_parser.add_argument(*verbose['args'], **verbose['kwargs'])
view3d_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])

#===============================================================================
# export
#===============================================================================
export_parser = subparsers.add_parser('export', description="Export segmentation as various file formats", help="export as file")
export_parser.add_argument('sff_file', help="any SFF file")
export_parser.add_argument(*primary_descriptor['args'], **primary_descriptor['kwargs'])
export_parser.add_argument(*normals_off['args'], **normals_off['kwargs'])
export_parser.add_argument(*mask_value['args'], **mask_value['kwargs'])
export_parser.add_argument(*verbose['args'], **verbose['kwargs'])
export_parser.add_argument(*output_path['args'], **output_path['kwargs'])

# get the full list of tools from the Parser object
tool_list = Parser._actions[1].choices.keys()
tool_list += ['parser', 'sffreader', 'meshreader', 'surfreader', 'modreader', 'omero_wrapper', 'meshtools', 'stlreader', 'mapreader']

# tests
test_help = "one or none of the following: {}".format(", ".join(tool_list))
tests_parser = subparsers.add_parser('tests', description="Run unit tests", help="run unit tests")
tests_parser.add_argument('tool', nargs='*', default='all', help=test_help)
tests_parser.add_argument('-v', '--verbosity', default=1, type=int, help="set verbosity; valid values: %s [default: 0]" % ", ".join(map(str, verbosity_range)))

test_parser = subparsers.add_parser('test', description="Run unit tests", help="run unit tests")
test_parser.add_argument('tool', nargs='*', default='all', help=test_help)
test_parser.add_argument('-v', '--verbosity', default=1, type=int, help="set verbosity; valid values: %s [default: 0]" % ", ".join(map(str, verbosity_range)))


def parse_args(_args):
    """
    Parse and check command-line arguments
    
    Subcommand handlers defined in __main__.py (e.g. handle_conver(...)) should not have to check arguments for consistency
    
    :param list _args: list of arguments
    :return: parsed arguments
    :rtype: `argparse.Namespace`
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
        if _args[0] in Parser._actions[1].choices.keys():
            exec('{}_parser.print_help()'.format(_args[0]))
            sys.exit(0)
    
    # parse args         
    args = Parser.parse_args(_args)
    
    # configs
#     configs = get_configs()
    
    # credentials
#     if args.subcommand in ['list', 'attachroi', 'delroi']:
#         if args.local 
    
    
    # check values
    # attachroi
#     if args.subcommand == 'attachroi':
#         # ensure that at least -x, -y, or -z is specified
#         if not args.x_image_id and not args.y_image_id and not args.z_image_id:
#             raise ValueError("attachroi requires at least one image ID (-x, -y, or -z)")

    # createroi
    if args.subcommand == 'createroi':
        # ensure valid primary_descriptor
        if args.primary_descriptor:
            try:
                assert args.primary_descriptor in ['threeDVolume', 'contourList', 'meshList', 'shapePrimitiveList']
            except AssertionError:
                print_date('Invalid value for primaryDescriptor: %s' % args.primary_descriptor)
                return None
        
        # ensure valid transparencey
        if args.transparency:
            try:
                assert 0 <= args.transparency <= 1
            except AssertionError:
                print_date("Invalid value for transparency: {}; should be between 0 and 1 (inclusive)".format(args.transparency))
                return None
        
        # ensure mask value is an integer (or long)
        if args.mask_value:
            try:
                assert isinstance(args.mask_value, int) or isinstance(args.mask_value, long)
            except AssertionError:
                print_date("Non-integer for mask value")
                return None 
    # delroi
    elif args.subcommand == 'delroi':
        # ensure that we have either an image or ROI ID
        if not args.image_id and not args.roi_id:
            raise ValueError('Missing either image (-i) or ROI (-r) ID')
        
        # ensure that both image and ROI ID are not set simultaneously
        if args.image_id and args.roi_id:
            raise ValueError('Only set one of image (-i) or ROI (-r) ID; not both')
    
    # tests
    elif args.subcommand == 'tests' or args.subcommand == 'test':        
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
                raise ValueError("Verbosity should be in %s-%s: %s given" % (verbosity_range[0], verbosity_range[-1], args.verbosity))
    
    # view3d
    elif args.subcommand == 'view3d':
        if args.primary_descriptor:
            try:
                assert args.primary_descriptor in ['threeDVolume', 'contourList', 'meshList', 'shapePrimitiveList']
            except:
                raise ValueError('Invalid value for primaryDescriptor: %s' % args.primary_descriptor)
        
        # ensure valid transparencey
        assert 0 <= args.transparency <= 1
        
        # ensure valid background colours
        assert 0 <= args.background_colour[0] <= 1
        assert 0 <= args.background_colour[1] <= 1
        assert 0 <= args.background_colour[2] <= 1
        
        # view contours
        # don't specify -A and -X, -Y, and/or -Z
        if args.all_contours or args.x_contours or args.y_contours or args.z_contours:
            assert (args.all_contours and not (args.x_contours or args.y_contours or args.z_contours)) or (not args.all_contours and (args.x_contours or args.y_contours or args.z_contours))
        
        # cube axes validity
        assert (0 <= args.cube_axes <= 4) or (args.cube_axes is None)
    
    return args 