#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sfftkplus.sffplus


"""

from __future__ import division

from getpass import getpass
import re
import os
import sys

import h5py

from .core.configs import get_configs
from sfftk.core.print_tools import print_date, print_static


# global dictionary containing persistent configs
configs = get_configs()


__author__  = "Paul K. Korir, PhD"
__email__   = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__    = "2017-02-17"


# auxilliary function
def get_credentials(args, configs):
    """Get credentials to use for OMERO connections from either args or configs
    
    :param args: command line argments
    :type args: ``argparse.Namespace``
    :param dict configs: dictionary of configurations read from text config file
    :return str user: username
    :return str password: password
    :return str host: hostname
    :return str port: port string
    """
    cw = configs['CONNECT_WITH']
    # user
    if 'OMERO_{}_USER'.format(cw) in configs:
        user = configs['OMERO_{}_USER'.format(cw)]
    elif args.user:
        user = args.user
    
    # password
    if 'OMERO_{}_PASSWORD'.format(cw) in configs:
        password = configs['OMERO_{}_PASSWORD'.format(cw)]
    elif args.password:
        password = args.password
    else:
        password = getpass()
   
    # host
    if 'OMERO_{}_HOST'.format(cw) in configs:
        host = configs['OMERO_{}_HOST'.format(cw)]
    elif args.host:
        host = args.host
       
    # port
    if 'OMERO_{}_PORT'.format(cw) in configs:
        port = configs['OMERO_{}_PORT'.format(cw)]
    elif args.port:
        port = args.port
    
    return user, password, host, port

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
    else: # check for image_ids in the roi file
        image_ids['x'] = roi_seg.header.front_id
        image_ids['y'] = roi_seg.header.right_id
        image_ids['z'] = roi_seg.header.top_id
    # make sure we have some image_ids to use
    if len(image_ids) == 0:
        raise ValueError("No image identifiers specified. Aborting...")
    return image_ids

def handle_list(args):
    from .omero.handlers import OMEROConnection
    with OMEROConnection(args) as connection:
        connection.list()
    return 0

def handle_attachroi(args):
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
    with OMEROConnection(args) as connection:
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
                continue # non-fatal
    return 0

def handle_delroi(args):
    from .omero.handlers import OMEROConnection
    with OMEROConnection(args) as connection:
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

def handle_createroi(args):
    """
    Handle `createroi` subcommand
    
    :param args: parsed arguments
    :type args: `argparse.Namespace`
    """
    import schema
    # convert an EMDB-SFF file to an ROI file
    if re.match(r'.*\.(sff|hff)$', args.sff_file, re.IGNORECASE):
        sff_seg = schema.SFFPSegmentation(args.sff_file)
    else:
        print_date("Unsupported file type: {}".format(args.sff_file))
        return 1  
    # convert segments to VTK meshes
    vtk_seg = sff_seg.as_vtk(args)
    # slice to get contours
    vtk_seg.slice()
    # convert to ROI using sfftkplus.schema.roi
    roi_seg = vtk_seg.as_roi()
    # export to file
    roi_seg.export(args.output)
    return 0

def handle_view3d(args):
    """
    Handle `view3d` subcommand
    
    :param args: parsed arguments
    :type args: `argparse.Namespace`
    """
    from sfftkplus import schema
    # convert an EMDB-SFF file to an ROI file
    if re.match(r'.*\.(sff|hff)$', args.sff_file, re.IGNORECASE):
        sff_seg = schema.SFFPSegmentation(args.sff_file)
    else:
        print_date("Unsupported file type: {}".format(args.sff_file))
        return 1
    # convert segments to VTK meshes
    vtk_seg = sff_seg.as_vtk(args)
    # render as 3D
    vtk_seg.render()
    return 0

def handle_export(args):
    if args.verbose:
        print_date("Converting segments in {} to VTP files".format(args.sff_file))
    from sfftkplus import schema
    if re.match(r'.*\.(sff|hff)$', args.sff_file, re.IGNORECASE):
        sff_seg = schema.SFFPSegmentation(args.sff_file)
    else:
        print_date("Unsupported file type: {}".format(args.sff_file))
        return 1
    vtk_seg = sff_seg.as_vtk(args)
    out_fn = os.path.basename(".".join(args.sff_file.split('.')[:-1]))
    vtk_seg.export(out_fn, args)
    return 0

def main():
    try:
        from sfftkplus.core.parser import parse_args
        args = parse_args(sys.argv[1:])
        # missing args
        if not args:
            return 1
        # subcommands
        if args.subcommand == 'createroi':
            return handle_createroi(args)
        elif args.subcommand == 'attachroi':
            return handle_attachroi(args)
        elif args.subcommand == 'delroi':
            return handle_delroi(args)
        elif args.subcommand == "view3d":
            return handle_view3d(args)
        elif args.subcommand == "list":
            return handle_list(args)
        elif args.subcommand == "export":
            return handle_export(args)
    except KeyboardInterrupt:
        # handle keyboard interrupt
        return 0

    return 0

if __name__ == "__main__":
    sys.exit(main())