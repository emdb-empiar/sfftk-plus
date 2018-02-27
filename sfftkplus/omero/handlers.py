# -*- coding: utf-8 -*-
"""
handler.py
============

Module OMERO handlers


TODO:
* proper implementation of __exit__ for OMEROConnectionHandler class


Version history:
0.0.1, 2016-02-22, First incarnation

"""

from __future__ import division

import re
import sys
import time

import omero.gateway
from omero.model import RoiI  # @UnresolvedImport
from omero.rtypes import rstring
import omero.sys
from sfftk.core.print_tools import print_date
from sfftkplus.omero import primitives


# from matplotlib import image
__author__ = "Paul K. Korir, PhD"
__email__ = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__ = "2016-02-22"

# view width
VIEW_WIDTH = 140


class OMEROROI(RoiI):
    """Wrapper around on OMERO ROI that inherits from ``omero.model.RoiI``"""
    def __init__(self, name, image, orientation=None, *args, **kwargs):
        """
        Initialiser of an ROI by name and image object
        
        :param str name: the name of this ROI
        :param `omero.gateway._ImageWrapper` image: the associated image
        """
        super(OMEROROI, self).__init__(*args, **kwargs)
        if isinstance(name, str):
            self.setName(name)
        if isinstance(image, omero.gateway._ImageWrapper):
            self.setImage(image)
            self.__image = image
        if orientation is not None:
            directions = ['x', 'y', 'z']
            try:
                assert orientation in directions
            except AssertionError:
                print_date("Invalid orientation: {}; should be in {}".format(orientation, ", ".join(directions)))
                sys.exit(1)
            self._orientation = orientation
        # try to determine the orientation from the image name
        else:
            if re.match(r'.*\-front\..*$', image.name, re.IGNORECASE):
                self._orientation = 'x'
            elif re.match(r'.*\-side\..*$', image.name, re.IGNORECASE):
                self._orientation = 'y'
            elif re.match(r'.*\-top\..*$', image.name, re.IGNORECASE):
                self._orientation = 'z'
            else:
                print_date("Unable to determine ROI orientation. Ensure image name has one of the form '*-front.*' (x), '*-side.*' (y) or '*-top.*' (z).")
                sys.exit(1)

    def getId(self):
        return super(OMEROROI, self).getId()

    def getName(self):
        return super(OMEROROI, self).getName()

    def setName(self, name):
        super(OMEROROI, self).setName(rstring(name))

    def getImage(self):
        return super(OMEROROI, self).getImage()

    def setImage(self, image):
        super(OMEROROI, self).setImage(image._obj)

    def getProject(self):
        return self.__image.getProject().getName()

    def getDatasets(self):
        return self.__image.listParents()

    @property
    def orientation(self):
        return self._orientation

    def addShape(self, shape):
        from omero.model import Shape  # @UnresolvedImport
        if isinstance(shape, Shape):
            super(OMEROROI, self).addShape(shape)
        else:
            raise TypeError('non-Shape argument: %s' % shape.__class__)

    def load_data(self, orientation_value, oriented_segments, segment_colours):
        for segment_id, contours in oriented_segments.iteritems():
            for contour in contours:
                theZ = orientation_value
                polygon = primitives.Polygon(theZ=theZ, theT=0)
                if self._orientation == 'x':
                    points = [(p.get_y(), p.get_z()) for p in contour.p]
                    polygon.setPoints(points, offsetYFrom=int(self.__image.getSizeX()))
                elif self._orientation == 'y':
                    points = [(p.get_x(), p.get_z()) for p in contour.p]
                    polygon.setPoints(points, offsetXFrom=int(self.__image.getSizeX()), offsetYFrom=int(self.__image.getSizeY()))
                elif self._orientation == 'z':
                    points = [(p.get_x(), p.get_y()) for p in contour.p]
                    polygon.setPoints(points, swapXY=True)
                polygon.setStrokeWidth(0.22, 'PIXEL')
                r, g, b, a = segment_colours[segment_id]
                polygon.setFillColor(r, g, b, A=a)
                polygon.setStrokeColor(0, 1, 0)
                polygon.setTextValue(segment_id)
                polygon.setFontSize(2, 'POINT')
                polygon.setFontStyle("Bold")
                self.addShape(polygon)


class OMEROROIList(object):
    def __init__(self, roi_seg, orientation, image, args):
        self._roi_seg = roi_seg
        self._orientation = orientation
        self._image = image
        self._args = args
        self._omero_rois = self._make_rois()

    def __iter__(self):
        return iter(self._omero_rois)

    def _make_rois(self):
        omero_rois = list()
        # we are in an orientation, say 'x'
        # 'x' has x-values as key and a value of a dictionary of segment_ids to contours
        # colour as rgba from colour 'key'
        segment_colours = self._roi_seg.segment_colours
        for orientation_value, oriented_segments in self._roi_seg.oriented_segments[self._orientation].iteritems():
            omero_roi = OMEROROI(name="{} = {}".format(self._orientation, orientation_value), image=self._image, orientation=self._orientation)
            omero_roi.load_data(orientation_value, oriented_segments, segment_colours)
            omero_rois.append(omero_roi)
        return omero_rois

    def __len__(self):
        return len(self._omero_rois)

    def __getitem__(self, index):
        return self._omero_rois[index]


class ImageView(object):
    width = VIEW_WIDTH
    def __init__(self, images):
        self._images = images

    def __iter__(self):
        return iter(self._images)

    def __str__(self):
        string = ""
        string += "*" * self.width + "\n"
        string += "*** AVAILABLE IMAGES *** ".center(self.width) + "\n"
        string += "*" * self.width + "\n"
        header = "" + \
            "ImageID".ljust(8) + \
            "ImageName".ljust(40) + \
            "Project".ljust(15) + \
            "Dataset".ljust(15) + \
            "X".rjust(6) + \
            "Y".rjust(6) + \
            "Z".rjust(6) + \
            u"pX (\U000003bcm)".rjust(10) + \
            u"pY (\U000003bcm)".rjust(10) + \
            u"pZ (\U000003bcm)".rjust(10) + \
            "#ROIs".rjust(6)
        string += header.encode('utf-8') + "\n"
        string += "-" * self.width + "\n"
        for image in self:
            string += "{:<8}{:<40}{:<15}{:<15}{:>6}{:>6}{:>6}{:>10.4}{:>10.4}{:>10.4}{:>6}\n".format(
                image.getId(),
                image.getName(),
                image.getProject().getName() if image.getProject() else "-",
                ", ".join(map(lambda x: x.getName(), image.listParents())),
                image.getSizeX(),
                image.getSizeY(),
                image.getSizeZ(),
                image.getPixelSizeX(),
                image.getPixelSizeY(),
                image.getPixelSizeZ(),
                image.getROICount(),
                )
        string += "-" * self.width + "\n"
        return string

    def __repr__(self):
        return str(self)


class ROIView(object):
    width = 140
    def __init__(self, rois):
        self._rois = rois

    def __iter__(self):
        return iter(self._rois)

    def __str__(self):
        string = ""
        string += "*" * self.width + "\n"
        string += "*** AVAILABLE ROIS ***".center(self.width) + "\n"
        string += "*" * self.width + "\n"
        header = "" + \
            "ROI ID".ljust(8) + \
            "ImageID".ljust(8) + \
            "ROI Name".ljust(20) + \
            "ImageName".ljust(40) + \
            "Project".ljust(15) + \
            "Dataset".ljust(15) + \
            "Primary shape".ljust(20)
        string += header.encode('utf-8') + "\n"
        string += "-" * self.width + "\n"
        for image, roi_list in self:
            for roi in roi_list:
                string += "{:<8}{:<8}{:<20}{:<40}{:<15}{:<15}{:<20}\n".format(
                    roi.getId().getValue(),
                    image.getId(),
                    roi.getName().getValue(),
                    image.getName(),
                    image.getProject().getName() if image.getProject() else "-",
                    ", ".join(map(lambda x: x.getName(), image.listParents())),
                    type(roi.getPrimaryShape()),
                    )
        string += "-" * self.width + "\n"
        return string

    def __repr__(self):
        return str(self)


class OMEROConnection(object):
    def __init__(self, args, configs):
        self.args = args
        self.configs = configs
        self.connect_with = self.configs['CONNECT_WITH']
        self.host = self.configs['OMERO_{}_HOST'.format(self.connect_with)]
        self.user = self.configs['OMERO_{}_USER'.format(self.connect_with)]
        self.password = self.configs['OMERO_{}_PASSWORD'.format(self.connect_with)]
        self.port = self.configs['OMERO_{}_PORT'.format(self.connect_with)]

    def __enter__(self):
        from omero.gateway import BlitzGateway  # @UnresolvedImport
        self.conn = BlitzGateway(
            host=self.host,
            username=self.user,
            passwd=self.password,
            port=self.port
            )
        self.connected = self.conn.connect()
        # failed to connect
        if not self.connected:
            print_date("Unable to connect to host '{}' on port {} using user '{}'".format(self.host, self.port, self.user))
            print_date("Check that server is up")
            sys.exit(1)
        # keepalive
        self.conn.c.enableKeepAlive(5)
        self.omero_user = self.conn.getUser()
        self.userId = self.omero_user.getId()
        self.updateService = self.conn.getUpdateService()
        self.roiService = self.conn.getRoiService()
        return self

    def __exit__(self, exc_type, exc_value, traceback):  # @UnusedVariable
        self.conn._closeSession()
        if self.args.verbose:
            print_date('Connection closed.')
        self.connected = False

    def list(self):
        """List images or ROIs
        
        :return int count: the number of images/ROIs
        """
        if self.args.images:
            print_date("Listing images...")
            images = self.images()
            if self.args.summary:
                pass
            else:
                print_date("Structuring output...")
                image_view = ImageView(images)
                print image_view
            return len(images)
        elif self.args.rois:
            print_date("Listing ROIs...")
            rois = self.rois()
            if self.args.summary:
                pass
            else:
                print_date("Structuring output...")
                roi_view = ROIView(rois)
                print roi_view
            return len(rois)

    @property
    def projects(self):
        if self.args.project is not None:
            projects = self.conn.searchObjects(["Project"], self.args.project)
            return projects
        else:
            return self.conn.listProjects(self.userId)

    def datasets(self, project=None):
        """List the datasets associated with the current user
        
        :param project: a project
        :param bool in_project: are these datasets contained within a project?
        """
        datasets = list()
        if self.args.dataset is not None:
            datasets += self.conn.searchObjects(["Dataset"], self.args.dataset)
        else:
            if project is not None:
                projects = self.conn.searchObjects(["Project"], project)
                for p in projects:
                    datasets += p.listChildren()
            else:
                params = omero.sys.ParametersI()  # @UndefinedVariable
                params.exp(self.userId)
                datasets += self.conn.getObjects("Dataset", params=params)
        return datasets

    def images(self, project=None, dataset=None):
        """The list of images associated with the current user
        
        If the image ID is specified only the required image is returned. (The project and dataset are ignored.)
        
        Otherwise:
        
        - If a project object is provided all images in all datasets in the project are returned. (The dataset is ignored.)
        - If a project object and dataset object are provided then only those images in the project and dataset are return.
        - If no project object is provided but a dataset object is provided then only those images in the dataset are returned. 
         
        :param str project: OMERO project name
        :param str dataset: OMERO dataset name
        :return list images: a list of OMERO ``Image`` objects
        """
        # assertions
        images = list()
        print_date("Retrieving images...")
        if self.args.image_id is not None:
            try:
                assert isinstance(self.args.image_id, int) or isinstance(self.args.image_id, long)
            except AssertionError:
                print_date("Invalid type for image ID: {}".format(type(self.args.image_id)))
                sys.exit(1)
            image = self.getImage(self.args.image_id)
            if image is not None:
                images.append(image)
        elif self.args.image_name is not None:
            images = self.conn.searchObjects(["Image"], self.args.image_name)
        else:
            if project is not None:  # project specified
                print_date("Searching for images in project '{}'".format(project))
                # get all projects matching
                projects = self.conn.searchObjects(["Project"], project)
                # get all datasets in projects matching
                datasets_in_projects = dict()
                for p in projects:
                    for d in p.listChildren():
                        datasets_in_projects[d.getName()] = d
                print_date("Found {} datasets in project '{}'".format(len(datasets_in_projects), project))
                # dataset specified
                if dataset is not None:
                    print_date("Searching for images in dataset '{}'".format(dataset))
                    if dataset in datasets_in_projects.keys():
                        images += datasets_in_projects[dataset].listChildren()
                else:  # dataset not specified
                    print_date("Searching for images in all {} datasets".format(len(datasets_in_projects)))
                    for dataset in datasets_in_projects.keys():
                        images += datasets_in_projects[dataset].listChildren()
            else:  # project not specified
                # dataset specified
                if dataset is not None:
                    print_date("Searching for images in dataset '{}'".format(dataset))
                    datasets = self.conn.searchObjects(["Dataset"], dataset)
                    for dataset in datasets:
                        images += dataset.listChildren()
                else:
                    datasets = self.datasets()
                    print_date("Searching for images in all {} datasets".format(len(datasets)))
                    for dataset in datasets:
                        images += dataset.listChildren()
        print_date("Found {} image(s).".format(len(images)))
        return images

    def getImage(self, image_id):
        """Get the image with the image ID specified on the command line
        
        :param int image_id: command line arguments
        :return image: an image
        :rtype image: ``OMEROImage``
        """
        return self.conn.getObject("Image", image_id)

    def rois(self, project=None, dataset=None):
        """Get an iterator over the ROIs associated with the specified image ID
        
        :param int image_id: image ID
        """
        rois = list()
        print_date("Retrieving ROIs...")
        if self.args.image_id is not None:
            return self.getROIs(self.args.image_id)
        else:
            for image in self.images(project, dataset):
                if image.getROICount() > 0:
                    rois.append((image, self.getROIs(image.getId())))
        roi_count = sum(map(lambda r: len(r[1]), rois))
        print_date("Found {:,} ROIs in {:,} images.".format(roi_count, len(rois)))
        return rois

    def getROIs(self, image_id):
        result = self.roiService.findByImage(image_id, None)
        return result.rois

    def attachRois(self, omero_rois):
        """Attach the rois from the iterable"""
        non_rois = filter(lambda r: not isinstance(r, OMEROROI), omero_rois)
        try:
            assert len(non_rois) == 0
        except AssertionError:
            print_date("Found {:,} non-ROI objects".format(len(non_rois)))
            return 1
        for roi in omero_rois:
            self.saveRoi(roi)
        return 0

    # save
    def saveRoi(self, roi):
        """Save the given ROI
        
        :param roi: an ROI object
        :type roi: `omero.model.Roi`
        """
        import Ice
        try:
            self.updateService.saveObject(roi)
        except Ice.MemoryLimitException as e:  # @UndefinedVariable
            print_date(str(e))
            sys.exit(1)

    # delete
    def deleteRoi(self, roi_id):
        """
        Delete the given ROI
        
        :param roi: an ROI object
        :type roi: `omero.model.Roi`
        """
        from omero.callbacks import CmdCallbackI  # @UnresolvedImport

        handle = self.conn.deleteObjects("Roi", [roi_id], deleteAnns=True, deleteChildren=True)
        callback = CmdCallbackI(self.conn.c, handle)

        while not callback.block(500):
            if self.args.verbose:
                print_date(".", newline=False, incl_date=False)
                time.sleep(2)

        callback.close(True)
