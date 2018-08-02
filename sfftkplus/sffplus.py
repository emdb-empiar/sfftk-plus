#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sfftkplus.sffplus


"""

from __future__ import division

import os
import re
import sys

from sfftk.core.print_tools import print_date, print_static
from sfftk.sff import _module_test_runner, _discover_test_runner

__author__ = "Paul K. Korir, PhD"
__email__ = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__ = "2017-02-17"
__updated__ = "2018-02-21"


# auxilliary function
# def get_credentials(args, configs):
#     """Get credentials to use for OMERO connections from either args or configs
#
#     :param args: command line argments
#     :type args: ``argparse.Namespace``
#     :param dict configs: dictionary of configurations read from text config file
#     :return str user: username
#     :return str password: password
#     :return str host: hostname
#     :return str port: port string
#     """
#     cw = configs['CONNECT_WITH']
#     # user
#     if 'OMERO_{}_USER'.format(cw) in configs:
#         user = configs['OMERO_{}_USER'.format(cw)]
#     elif args.user:
#         user = args.user
#
#     # password
#     if 'OMERO_{}_PASSWORD'.format(cw) in configs:
#         password = configs['OMERO_{}_PASSWORD'.format(cw)]
#     elif args.password:
#         password = args.password
#     else:
#         password = getpass()
#
#     # host
#     if 'OMERO_{}_HOST'.format(cw) in configs:
#         host = configs['OMERO_{}_HOST'.format(cw)]
#     elif args.host:
#         host = args.host
#
#     # port
#     if 'OMERO_{}_PORT'.format(cw) in configs:
#         port = configs['OMERO_{}_PORT'.format(cw)]
#     elif args.port:
#         port = args.port
#
#     return user, password, host, port

def get_image_ids(roi_seg, args):
    """Get the image ids to which ROIs will be attached
    
    :param roi_seg: an ROI segmentation object
    :type roi_seg: ``sfftkplus.formats.roi.ROISegmentation``
    :param args: command line arguments
    :type args: ``argparse.Namespace``
    :return dict image_ids: dictionary of image ids index by orientation ('x', 'y', 'z')  
    """
    image_ids = dict()
    if args.x_image_id or args.y_image_id or args.z_image_id:
        if args.x_image_id:
            image_ids['x'] = args.x_image_id
        if args.y_image_id:
            image_ids['y'] = args.y_image_id
        if args.z_image_id:
            image_ids['z'] = args.z_image_id
    else:  # check for image_ids in the roi file
        image_ids['x'] = roi_seg.header.front_id
        image_ids['y'] = roi_seg.header.right_id
        image_ids['z'] = roi_seg.header.top_id
    #  make sure we have some image_ids to use
    if len(image_ids) == 0:
        raise ValueError("No image identifiers specified. Aborting...")
    return image_ids


def handle_list(args, configs):
    """
    Handle `list` subcommand
    
    :param args: parsed arguments
    :type args: ``argparse.Namespace``
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int status: status
    """
    from .omero.handlers import OMEROConnection
    with OMEROConnection(args, configs) as connection:
        connection.list()
    return 0


def handle_attachroi(args, configs):
    """
    Handle `attachroi` subcommand
    
    :param args: parsed arguments
    :type args: ``argparse.Namespace``
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int status: status
    """
    from .formats.roi import ROISegmentation
    from .omero.handlers import OMEROConnection

    if re.match(r'.*\.roi$', args.roi_file, re.IGNORECASE):
        if args.verbose:
            print_date("Reading ROIs from %s..." % args.roi_file)
        roi_seg = ROISegmentation(args.roi_file)
    else:
        print_date("Unkown file type '%s'. Should be valid ROI file. Aborting..." % args.roi_file)
        return 1
    # image_ids
    image_ids = get_image_ids(roi_seg, args)
    # open the connection
    # iterate over orientations
    #    get image id for this orientation
    #    check image exists
    #    create and load roi object
    #    attach rois
    with OMEROConnection(args, configs) as connection:
        for orientation in roi_seg.oriented_segments:
            # get the image_id for this orientation
            image_id = image_ids[orientation]
            if args.verbose:
                print_date("Checking whether %s image of ID %s exists..." % (orientation, image_id), newline=False)
            # get the image to which we will write ROIs
            image = connection.getImage(image_id)
            if image:
                if args.verbose:
                    print_date("OK", incl_date=False)
                # convert the oriented segments to rois
                omero_rois = roi_seg.as_omero_rois(orientation, image, args)
                # load the rois to OMERO
                if args.verbose:
                    print_date("Attaching ROIs...", newline=False)
                status = connection.attachRois(omero_rois)
                if status == 0:
                    if args.verbose:
                        print_date("OK", incl_date=False)
                else:
                    if args.verbose:
                        print_date("FAIL", incl_date=False)
                    return status
            else:
                if args.verbose:
                    print_date("FAIL", incl_date=False)
                continue  # non-fatal
    return 0


def handle_delroi(args, configs):
    """
    Handle `delroi` subcommand
    
    :param args: parsed arguments
    :type args: ``argparse.Namespace``
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int status: status
    """
    from .omero.handlers import OMEROConnection
    with OMEROConnection(args, configs) as connection:
        if args.roi_id:
            if args.verbose:
                print_date("Deleting ROI %s" % args.roi_id)
            connection.deleteRoi(args.roi_id)
        elif args.image_id:
            rois = connection.rois(args.image_id)
            roi_count = len(rois)
            for roi in rois:
                if args.verbose:
                    print_static("Deleting ROI %s" % roi.getId().getValue())
                connection.deleteRoi(roi.getId().getValue())
            print_static("\n", incl_date=False)
            print_date("Deleted {} ROIs".format(roi_count))
        else:
            print_date("Please specify an ROI ID. Search using 'sffp list --rois [--image-id <image_id>]'")
    return 0


def handle_createroi(args, configs):
    """
    Handle `createroi` subcommand
    
    :param args: parsed arguments
    :type args: ``argparse.Namespace``
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int exit_status: exit status
    """
    # convert an EMDB-SFF file to an ROI file
    if re.match(r'.*\.(sff|hff)$', args.sff_file, re.IGNORECASE):
        from .schema import SFFPSegmentation
        if args.verbose:
            print_date("Reading in EMDB-SFF file {}".format(args.sff_file))
        sff_seg = SFFPSegmentation(args.sff_file)
        # convert segments to VTK meshes
        if args.verbose:
            print_date("Converting EMDB-SFF segments to VTK meshes")
        vtk_seg = sff_seg.as_vtk(args, configs)
        # slice to get contours
        if args.verbose:
            print_date("Slicing segmentation to get ROI contours...")
        vtk_seg.slice()
        # convert to ROI using sfftkplus.schema.roi
        if args.verbose:
            print_date("Converting to ROI using roi.xsd...")
        roi_seg = vtk_seg.as_roi(args, configs)
        # export to file
        if args.verbose:
            print_date("Writing output to {}".format(args.output))
        exit_status = roi_seg.export(args.output, args)
        if args.verbose:
            print_date("Done")
    elif re.match(r'.*\.roi$', args.sff_file, re.IGNORECASE):
        from .formats import roi
        if args.verbose:
            print_date("Reading in ROI file {}".format(args.sff_file))
        roi_seg = roi.ROISegmentation(args.sff_file)
        if args.reset_ids:
            if args.verbose:
                print_date("Resetting IDs...")
            roi_seg.header.reset_ids(args, configs)
        # export to file
        if args.verbose:
            print_date("Writing output to {}".format(args.output))
        exit_status = roi_seg.export(args.output)
        if args.verbose:
            print_date("Done")
    else:
        print_date("Unsupported file type: {}".format(args.sff_file))
        exit_status = os.EX_DATAERR
    return exit_status


def handle_view3d(args, configs):
    """
    Handle `view3d` subcommand
    
    :param args: parsed arguments
    :type args: ``argparse.Namespace``
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int status: status
    """
    from sfftkplus import schema
    # convert an EMDB-SFF file to an ROI file
    if re.match(r'.*\.(sff|hff)$', args.sff_file, re.IGNORECASE):
        sff_seg = schema.SFFPSegmentation(args.sff_file)
    else:
        print_date("Unsupported file type: {}".format(args.sff_file))
        return 1
    # convert segments to VTK meshes
    vtk_seg = sff_seg.as_vtk(args, configs)
    # render as 3D
    vtk_seg.render()
    return 0


def handle_export(args, configs):
    """
    Handle `export` subcommand
    
    :param args: parsed arguments
    :type args: ``argparse.Namespace``
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int status: status
    """
    if args.verbose:
        print_date("Converting segments in {} to VTP files".format(args.sff_file))
    from . import schema
    if re.match(r'.*\.(sff|hff)$', args.sff_file, re.IGNORECASE):
        sff_seg = schema.SFFPSegmentation(args.sff_file)
    else:
        print_date("Unsupported file type: {}".format(args.sff_file))
        return 1
    vtk_seg = sff_seg.as_vtk(args, configs)
    out_fn = os.path.basename(".".join(args.sff_file.split('.')[:-1]))
    vtk_seg.export(out_fn, args)
    return 0


def handle_configs(args, configs):
    """
    Handle `configs` subcommand
    
    :param args: parsed arguments
    :type args: ``argparse.Namespace``
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int status: status
    """
    if args.config_subcommand == "list":
        from sfftk.core.configs import list_configs
        return list_configs(args, configs)
    elif args.config_subcommand == "get":
        from sfftk.core.configs import get_configs
        return get_configs(args, configs)
    elif args.config_subcommand == "set":
        from sfftk.core.configs import set_configs
        return set_configs(args, configs)
    elif args.config_subcommand == "del":
        from sfftk.core.configs import del_configs
        return del_configs(args, configs)
    elif args.config_subcommand == "clear":
        from sfftk.core.configs import clear_configs
        return clear_configs(args, configs)
    return 0


def handle_tests(args, configs):
    """Handle `test` subcommand
    
    :param args: parsed arguments
    :type args: `argparse.Namespace`
    :param configs: configurations object
    :type config: ``sfftk.core.configs.Congif``
    :return int status: status
    """
    if isinstance(args.tool, str):
        #         from .unittests import test_main
        #         _module_test_runner(test_main, args)
        _discover_test_runner("sfftkplus.unittests", args)
    else:
        #         if 'main' in args.tool:
        #             from .unittests import test_main
        #             _module_test_runner(test_main, args)
        if 'core' in args.tool:
            from .unittests import test_core
            _module_test_runner(test_core, args)
        if 'formats' in args.tool:
            from .unittests import test_formats
            _module_test_runner(test_formats, args)
        if 'schema' in args.tool:
            from .unittests import test_schema
            _module_test_runner(test_schema, args)
        if 'readers' in args.tool:
            from .unittests import test_readers
            _module_test_runner(test_readers, args)
        if 'omero' in args.tool:
            from .unittests import test_omero
            _module_test_runner(test_omero, args)
    return 0


def main():
    try:
        from sfftkplus.core.parser import parse_args
        args, configs = parse_args(sys.argv[1:])
        # missing args
        if not args:
            return 1
        # subcommands
        if args.subcommand == "config":
            return handle_configs(args, configs)
        elif args.subcommand == 'createroi':
            return handle_createroi(args, configs)
        elif args.subcommand == 'attachroi':
            return handle_attachroi(args, configs)
        elif args.subcommand == 'delroi':
            return handle_delroi(args, configs)
        elif args.subcommand == "view3d":
            return handle_view3d(args, configs)
        elif args.subcommand == "list":
            return handle_list(args, configs)
        elif args.subcommand == "export":
            return handle_export(args, configs)
        elif args.subcommand == "tests":
            return handle_tests(args, configs)
    except KeyboardInterrupt:
        # handle keyboard interrupt
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
