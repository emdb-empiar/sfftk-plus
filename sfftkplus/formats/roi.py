#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sfftkplus.formats.roi



Copyright 2017 EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
You may obtain a copy of the License at 

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software 
distributed under the License is distributed on an 
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, 
either express or implied. 

See the License for the specific language governing permissions 
and limitations under the License.
"""

from __future__ import division

import os
import sys

import psycopg2

from sfftk.core.print_tools import print_date
from sfftk.core.utils import rgba_to_hex
from ..formats.base import Segmentation, Header, Annotation, Segment, Contours
from ..readers import roireader
from ..schema import roi

__author__ = "Paul K. Korir, PhD"
__email__ = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__ = "2017-04-11"

ORIENTATIONS = ['x', 'y', 'z']


def get_image_id(cursor, image_name_root, view, ext='map', quick_pick=None):
    """Obtain the image IDs for top, front and right images by EMDB accession code
    
    :param cursor: cursor to postgres connection
    :type cursor: psycopg2.Cursor
    :param str image_name_root: accession code in lowercase e.g. 'emd_1080'
    :param str view: the view of the image; either 'top', 'front' or 'right'
    :param str ext: extension of image file e.g. 'map'
    :return int image_id: the image ID or 0 for fail (no or multiple image IDs found)
    """
    try:
        assert isinstance(cursor, psycopg2.extensions.cursor)
    except AssertionError:
        print_date("Not pyscopg2.extensions.cursor object: {}".format(cursor))
        sys.exit(1)
    views = ['top', 'front', 'side']
    try:
        assert view in views
    except AssertionError:
        print_date("Invalid view: {}; should be one of: {}".format(view, ", ".join(views)))
        sys.exit(1)
    exts = ['map', 'mrc', 'rec']
    try:
        assert ext in exts  # supported file extensions
    except AssertionError:
        print_date("Invalid extension: {}; should be one of {}".format(ext, ", ".join(exts)))
        sys.exit(1)
    query_string = "select id from image where image.name like '{}-{}.%'".format(image_name_root, view)
    cursor.execute(query_string)
    rows = cursor.fetchall()
    if rows:
        if len(rows) == 1:
            return rows[0][0]
        else:
            print_date("Multiple image IDs for {}-{}: {}".format(image_name_root, view, rows))
            if quick_pick is not None:
                print_date("Quick picking an ID from index {}".format(quick_pick))
                return rows[quick_pick][0]
            else:
                return 0
    else:
        print_date("No image IDs found for view '{}'".format(view))
        return 0


class ROIContours(Contours):
    def __init__(self, contours=None):
        if contours:
            self._contours = contours
            self._x_contours = self._contours.xContours.contour
            self._y_contours = self._contours.yContours.contour
            self._z_contours = self._contours.zContours.contour

    @classmethod
    def from_vtk(cls, vtk_seg):
        obj = cls()
        obj.x_contours = vtk_seg.x_contours
        obj.y_contours = vtk_seg.y_contours
        obj.z_contours = vtk_seg.z_contours
        return obj

    @property
    def x_contours(self):
        return self._x_contours

    @x_contours.setter
    def x_contours(self, value):
        assert isinstance(value, list)
        self._x_contours = value

    @property
    def y_contours(self):
        return self._y_contours

    @y_contours.setter
    def y_contours(self, value):
        assert isinstance(value, list)
        self._y_contours = value

    @property
    def z_contours(self):
        return self._z_contours

    @z_contours.setter
    def z_contours(self, value):
        assert isinstance(value, list)
        self._z_contours = value

    def convert(self, args, configs):
        # x contours
        xContours = roi.orientedContourType()
        # print(self.x_contours)
        for contour in self.x_contours:  # for each contour
            K = roi.contourType()
            for point_id, point in contour.iteritems():  # for each point in the contour
                p = roi.pointType()
                x, y, z = point
                p.set_id(point_id)
                p.set_x(x)
                p.set_y(y)
                p.set_z(z)
                K.add_p(p)
            xContours.add_contour(K)
        # y contours
        yContours = roi.orientedContourType()
        for contour in self.y_contours:
            K = roi.contourType()
            for point_id, point in contour.iteritems():
                p = roi.pointType()
                x, y, z = point
                p.set_id(point_id)
                p.set_x(x)
                p.set_y(y)
                p.set_z(z)
                K.add_p(p)
            yContours.add_contour(K)
        # z contours
        zContours = roi.orientedContourType()
        for contour in self.z_contours:
            K = roi.contourType()
            for point_id, point in contour.iteritems():
                p = roi.pointType()
                x, y, z = point
                p.set_id(point_id)
                p.set_x(x)
                p.set_y(y)
                p.set_z(z)
                K.add_p(p)
            zContours.add_contour(K)
        return xContours, yContours, zContours


class ROISegment(Segment):
    def __init__(self, segment=None):
        if segment:
            self._segment = segment
            self._id = segment.id
            self._colour = segment.colour.red, segment.colour.green, segment.colour.blue, segment.colour.alpha
            self._contours = ROIContours(segment)
            self._oriented_contours = self._compute_oriented_contours()

    @classmethod
    def from_vtk(cls, vtk_seg):
        obj = cls()
        obj.id = vtk_seg.id
        obj.colour = vtk_seg.colour
        obj.contours = ROIContours.from_vtk(vtk_seg.contours[0])
        return obj

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        assert id >= 0 and isinstance(value, int)
        self._id = value

    @property
    def colour(self):
        return self._colour

    @colour.setter
    def colour(self, value):
        # assertions?
        self._colour = value

    @property
    def contours(self):
        return self._contours

    @contours.setter
    def contours(self, value):
        assert isinstance(value, ROIContours)
        self._contours = value

    def _compute_oriented_contours(self, *args, **kwargs):
        oriented_contours = dict()
        oriented_contours['x'] = dict()
        oriented_contours['y'] = dict()
        oriented_contours['z'] = dict()
        # xContours
        for xContour in self.contours.x_contours:
            x = int(xContour.p[0].get_x())
            if x not in oriented_contours['x']:
                oriented_contours['x'][x] = [xContour]
            else:
                oriented_contours['x'][x] += [xContour]
        #  yContours
        for yContour in self.contours.y_contours:
            y = int(yContour.p[0].get_y())
            if y not in oriented_contours['y']:
                oriented_contours['y'][y] = [yContour]
            else:
                oriented_contours['y'][y] += [yContour]
        # zContours
        for zContour in self.contours.z_contours:
            z = int(zContour.p[0].get_z())
            if z not in oriented_contours['z']:
                oriented_contours['z'][z] = [zContour]
            else:
                oriented_contours['z'][z] += [zContour]
        return oriented_contours

    @property
    def oriented_contours(self):
        return self._oriented_contours

    def convert(self, args, configs):
        segment = roi.segmentType()
        segment.id = self.id
        segment.colour = roi.rgbaType()
        segment.colour.red, segment.colour.green, segment.colour.blue, segment.colour.alpha = self.colour
        xContours, yContours, zContours = self.contours.convert(args, configs)
        segment.set_xContours(xContours)
        segment.set_yContours(yContours)
        segment.set_zContours(zContours)
        return segment


class ROIAnnotation(Annotation):
    pass


class ROIHeader(Header):
    def __init__(self, segmentation=None):
        if segmentation:
            self._segmentation = segmentation
            if segmentation.image_ids:
                self._top_id = segmentation.image_ids.top
                self._front_id = segmentation.image_ids.front
                self._right_id = segmentation.image_ids.right
            else:
                self._top_id = None
                self._front_id = None
                self._right_id = None

    def reset_ids(self, args, configs, *args_, **kwargs_):
        try:
            if args.image_name_root is not None:
                cw = configs['CONNECT_WITH']  # either LOCAL or REMOTE
                # server settings
                conn_str = "dbname='{}' user='{}' password='{}' host='{}' port='{}'".format(
                    configs['IMAGE_DB_{}_NAME'.format(cw)],
                    configs['IMAGE_DB_{}_USER'.format(cw)],
                    configs['IMAGE_DB_{}_PASS'.format(cw)],
                    configs['IMAGE_DB_{}_HOST'.format(cw)],
                    configs['IMAGE_DB_{}_PORT'.format(cw)],
                )
                conn = psycopg2.connect(conn_str)
                cur = conn.cursor()
                self.top_id = get_image_id(cur, args.image_name_root, 'top', quick_pick=args.quick_pick)
                self.front_id = get_image_id(cur, args.image_name_root, 'front', quick_pick=args.quick_pick)
                self.right_id = get_image_id(cur, args.image_name_root, 'side', quick_pick=args.quick_pick)
                # sanity check
                assert self.top_id != self.front_id and self.right_id != self.front_id
            elif args.top_front_right is not None:
                self.top_id, self.front_id, self.right_id = args.top_front_right
            else:
                print_date(
                    "Neither -I/--image-name-root nor --top-front-right arguments not set. Image IDs will be excluded.")
                self.top_id = None
                self.front_id = None
                self.right_id = None
        except AssertionError:
            print_date("Invalid image IDs or image IDs not found. Did you use -I/--image-name-root option?")
            self.top_id = None
            self.front_id = None
            self.right_id = None

    @classmethod
    def from_vtk(cls, vtk_args, configs, *args, **kwargs):
        obj = cls()
        obj.reset_ids(vtk_args, configs, *args, **kwargs)
        return obj

    @property
    def top_id(self):
        return self._top_id

    @top_id.setter
    def top_id(self, value):
        try:
            assert (value > 0 and (isinstance(value, int) or isinstance(value, long)) or value is None)
        except AssertionError:
            print_date("Invalid value: {}".format(value))
        #             sys.exit(1)
        self._top_id = value

    @property
    def front_id(self):
        return self._front_id

    @front_id.setter
    def front_id(self, value):
        try:
            assert (value > 0 and (isinstance(value, int) or isinstance(value, long)) or value is None)
        except AssertionError:
            print_date("Invalid value: {}".format(value))
        #             sys.exit(1)
        self._front_id = value

    @property
    def right_id(self):
        return self._right_id

    @right_id.setter
    def right_id(self, value):
        try:
            assert (value > 0 and (isinstance(value, int) or isinstance(value, long)) or value is None)
        except AssertionError:
            print_date("Invalid value: {}".format(value))
        #             sys.exit(1)
        self._right_id = value

    def convert(self, *args, **kwargs):
        if all([self.top_id, self.front_id, self.right_id]):
            image_ids = roi.image_idsType()
            image_ids.set_top(self.top_id)
            image_ids.set_front(self.front_id)
            image_ids.set_right(self.right_id)
            return image_ids
        else:
            return None


class ROISegmentation(Segmentation):
    def __init__(self, fn=None, *args, **kwargs):
        self._fn = fn
        if fn:
            self.roi_seg = roireader.get_data(fn, *args, **kwargs)
            self._header = ROIHeader(self.roi_seg)
            self._segments = map(ROISegment, self.roi_seg.segment)
            self._oriented_segments, self._segment_colours = self._compute_oriented_segments()

    @classmethod
    def from_vtk(cls, vtk_seg, args, configs):
        obj = cls()
        obj.header = ROIHeader.from_vtk(args, configs)
        obj.segments = map(ROISegment.from_vtk, vtk_seg.segments)
        obj.convert(args, configs)
        return obj

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, value):
        assert isinstance(value, ROIHeader)
        self._header = value

    @property
    def segments(self):
        return self._segments

    @segments.setter
    def segments(self, value):
        assert isinstance(value, list)
        self._segments = value

    def convert(self, args, configs):
        self.roi_seg = roi.ROI()
        self.roi_seg.set_image_ids(self.header.convert())
        for segment in self.segments:
            self.roi_seg.add_segment(segment.convert(args, configs))

    def _compute_oriented_segments(self, *args, **kwargs):
        oriented_segments = {
            'x': dict(),
            'y': dict(),
            'z': dict(),
        }
        segment_colours = dict()
        for segment in self.segments:
            segment_colours[segment.id] = segment.colour
            for o in segment.oriented_contours:
                for ovalue in segment.oriented_contours[o]:
                    if ovalue not in oriented_segments[o]:
                        oriented_segments[o][ovalue] = dict()
                    if segment.id not in oriented_segments[o][ovalue]:
                        oriented_segments[o][ovalue][segment.id] = segment.oriented_contours[o][ovalue]
                    else:
                        oriented_segments[o][ovalue][segment.id] += segment.oriented_contours[o][ovalue]
        return oriented_segments, segment_colours

    @property
    def oriented_segments(self):
        return self._oriented_segments

    @property
    def segment_colours(self):
        return self._segment_colours

    def as_omero_rois(self, orientation, image, args):
        """Convert an ROISegmentation object to a set of OMERO ROIs"""
        from ..omero.handlers import OMEROROIList
        if args.verbose:
            print_date("Creating iterator of OMERO ROIs for %s..." % orientation, newline=False)
        omero_rois = OMEROROIList(self, orientation, image, args)
        if args.verbose:
            print_date("OK", incl_date=False)
        return omero_rois

    def _export_rois_json(self, fn, orientation, args, fill_alpha=1.0, stroke_alpha=1.0, font_size=2.0,
                          stroke_colour=(0, 1, 0), stroke_width=0.22):
        """Export ROIs for this orientation as JSON (instead of as XML)

        :param fn: the output file name root; the image ID (or orientation) and the slice value will be included
        :param orientation: character specifying the orientation; either 'x', 'y', or 'z'
        :param fill_alpha: the alpha value for the shape fill; default 1.0
        :param stroke_alpha: the alpha value for the shape stroke; default 1.0
        :param font_size: the font size for text
        :param stroke_colour: the colour of the shape stroke; default (0, 1, 0) (green)
        :param stroke_width: the width of the shape stroke ; default 0.22
        :return exit_status: the exit status (see Python's ``os`` module for details)
        """
        try:
            assert orientation in ORIENTATIONS
        except AssertionError as a:
            print_date("Invalid value for 'orientation'; should be in [{}]".format(', '.join(ORIENTATIONS)))
            print(str(a))
            return os.EX_DATAERR
        # ensure valid shape params
        try:
            assert 0 <= fill_alpha <= 1
        except AssertionError as a:
            print_date("Invalid value for 'fill_alpha' ({}); should be in [0-1]".format(fill_alpha))
            print(str(a))
            return os.EX_DATAERR
        try:
            assert 0 <= stroke_alpha <= 1
        except AssertionError as a:
            print_date("Invalid value for 'stroke_alpha' ({}); should be in [0-1]".format(stroke_alpha))
            print(str(a))
            return os.EX_DATAERR
        try:
            assert 0 <= font_size <= 30  # arbitrary bounds
        except AssertionError as a:
            print_date("Invalid value for 'font_size' ({}); should be in [0-30]".format(font_size))
            print(str(a))
            return os.EX_DATAERR
        try:
            assert 3 <= len(stroke_colour) <= 4
            if len(stroke_colour) == 3:
                r, g, b = stroke_colour
                a = 1
            elif len(stroke_colour) == 4:
                r, g, b, a = stroke_colour
            assert 0 <= r <= 1 and 0 <= g <= 1 and 0 <= b <= 1 and 0 <= a <= 1
        except AssertionError as a:
            print_date("Invalid value for 'stroke_colour' ({})".format(str(stroke_colour)))
            print(str(a))
            return os.EX_DATAERR
        try:
            assert 0 <= stroke_width <= 10
        except AssertionError as a:
            print_date("Invalid value for 'stroke_width' (); should be in [0-10]".format(stroke_width))
            print(str(a))
            return os.EX_DATAERR
        # stroke colour
        stroke_colour_ = rgba_to_hex(stroke_colour)
        import json
        # x contours
        grouped_contours = dict()
        # we want to group contours by slice value (o)
        # each json file will have all the contours for all the segments batched together
        # we need to store the colour for each segment and the contours for each segment
        for segment in self.segments:
            for contour in getattr(segment.contours, '{}_contours'.format(orientation)):
                # get the value for this orientation slice level
                if orientation == 'x':
                    o = int(contour[0][0])
                elif orientation == 'y':
                    o = int(contour[0][1])
                elif orientation == 'z':
                    o = int(contour[0][2])
                if o not in grouped_contours:
                    grouped_contours[o] = dict()
                if segment.id not in grouped_contours[o]:
                    grouped_contours[o][segment.id] = dict()
                    grouped_contours[o][segment.id]['colour'] = segment.colour
                    grouped_contours[o][segment.id]['contours'] = [contour]
                else:
                    grouped_contours[o][segment.id]['contours'] += [contour]
        # we write a json file for each slice
        for o, segment in grouped_contours.iteritems():
            shapes = list()
            for segment_id, segment_contours in segment.iteritems():
                colour = segment_contours['colour']
                for contour in segment_contours['contours']:
                    point_str = ''
                    for point_id, point in contour.iteritems():
                        if point_id == 0:
                            if orientation == 'x':
                                point_str += 'M {:.2f} {:.2f} '.format(point[1], point[2])
                            elif orientation == 'y':
                                point_str += 'M {:.2f} {:.2f} '.format(point[0], point[2])
                            elif orientation == 'z':
                                point_str += 'M {:.2f} {:.2f} '.format(point[1], point[0])
                        else:
                            if orientation == 'x':
                                point_str += 'L {:.2f} {:.2f} '.format(point[1], point[2])
                            elif orientation == 'y':
                                point_str += 'L {:.2f} {:.2f} '.format(point[0], point[2])
                            elif orientation == 'z':
                                point_str += 'L {:.2f} {:.2f} '.format(point[1], point[0])
                    # if the last contour point is the same as the first the it is closed
                    # if contour[point_id] == contour[0]:
                    point_str += 'z'
                    shapes.append({
                        "fontStyle": "Bold",
                        "fillAlpha": fill_alpha,
                        "strokeAlpha": stroke_alpha,
                        "id": None,
                        "points": point_str,
                        "fontSize": font_size,
                        "theZ": o,
                        "strokeColor": stroke_colour_,
                        "theT": 0,
                        "type": "Polygon",
                        "textValue": str(segment_id),
                        "strokeWidth": stroke_width,
                        "fillColor": rgba_to_hex(colour),  # if contour[point_id] == contour[0] else None,
                    })
            # write the shapes for this slice
            if orientation == 'x':
                if self.roi_seg.image_ids.front is not None:
                    ofn = fn.format(self.roi_seg.image_ids.front, o)
                else:
                    ofn = fn.format(orientation, o)
            elif orientation == 'y':
                if self.roi_seg.image_ids.right is not None:
                    ofn = fn.format(self.roi_seg.image_ids.right, o)
                else:
                    ofn = fn.format(orientation, o)
            elif orientation == 'z':
                if self.roi_seg.image_ids.top is not None:
                    ofn = fn.format(self.roi_seg.image_ids.top, o)
                else:
                    ofn = fn.format(orientation, o)
            # write out the JSON for this orientation and this slice
            with open(ofn, 'w') as f:
                json.dump([{"shapes": shapes}], f)
        return

    def _export_json(self, fn, args, *_args, **_kwargs):
        """Export ROIs as JSON formatted the way OMERO structures its JSONs

        :param fn: the root of the output file names; should contain two placeholders for the image ID (if it does not
        exist the orientation ('x', 'y', 'z') will be used instead; and the slice level value
        :param _args: positional arguments to be passed to the underlying ``_export_orientation_json()`` private method
        :param _kwargs: keyword arguments to be passed to the underlying ``_export_orientation_json() private method
        :return exit_status: the exit status (see Python's ``os`` module for details
        """
        for orientation in ORIENTATIONS:
            exit_status = self._export_rois_json(fn, orientation, args, *_args, **_kwargs)
        return exit_status

    def export(self, fn, args, *_args, **_kwargs):
        """Export ROIs as a file

        The file extension determines the file format:

        * ``.roi`` outputs an XML file according to ``roi.xsd``, which ships with ``sfftk-plus``

        * ``.json`` outputs one JSON file for each orientation slice. For example, if there are 20
        slices in the x-direction, 30 slices in the y-direction and 40 slices in the z-direction then
        there will be 90 JSON files in all

        :param fn: the output file name; the extension is important
        :param args: positional arguments to be passed on
        :param _kwargs: keyword arguments to be passed on
        :return: None
        """
        import re
        # create the path if it doesn't exist
        path = os.path.dirname(fn)
        if not os.path.exists(path):
            if args.verbose:
                print_date("Path not found: {}. It will be created".format(path))
            os.makedirs(path)
        else:
            if args.verbose:
                print_date("Path found: {}".format(path))
        if re.match(r".*\.roi$", fn, re.IGNORECASE):
            with open(fn, 'w') as f:
                self.roi_seg.set_image_ids(self.header.convert())
                version = _kwargs.get('version') if 'version' in _kwargs else "1.0"
                encoding = _kwargs.get('encoding') if 'encoding' in _kwargs else "UTF-8"
                f.write('<?xml version="{}" encoding="{}"?>\n'.format(version, encoding))
                self.roi_seg.export(f, 0)
            exit_status = os.EX_OK
        elif re.match(r".*\.json$", fn, re.IGNORECASE):
            # in the case of outputing JSON the provided filename is not the actual filename into which data will be
            # written; rather, the filename conveys: i) the filename base; ii) the output format
            # the fn_root var is constructed from the fn argument
            fn_root = '.'.join(fn.split('.')[:-1]) + '-{}-{}.json'
            return self._export_json(fn_root, args, *_args, **_kwargs)
