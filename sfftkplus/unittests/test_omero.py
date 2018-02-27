#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_omero_wrapper.py

Unit tests for OMERO handlers
'''

from __future__ import division

import os, sys
import shlex
import unittest

from omero.gateway import _ImageWrapper
from omero.model import RectangleI  # @UnresolvedImport

from . import TEST_DATA_PATH, _random_float, _random_integer
from ..core.parser import parse_args
from ..omero import handlers, primitives
from ..sffplus import get_image_ids


__author__ = 'Paul K. Korir, PhD'
__email__ = 'pkorir@ebi.ac.uk, paul.korir@gmail.com'
__date__ = '2016-06-13'


# utility functions
def colorInt(r, g, b, a):
    '''Function that return a signed 32-bit integer from RGBA tuple.
    
    0 <= r,g,b,a <= 1
    '''
    assert 0 <= r <= 1
    assert 0 <= g <= 1
    assert 0 <= b <= 1
    assert 0 <= a <= 1
    # the actual format is 'ARGB'!!!
    rgba_int = (int(r * 255) << 24) + (int(g * 255) << 16) + (int(b * 255) << 8) + int(a * 255)
    if rgba_int > 2147483647:
        rgba_int = -2147483648 + rgba_int % 2147483648
    return rgba_int


# handlers
class TestHandlers_OMEROConnection(unittest.TestCase):
    '''Test on using the opened connection'''
    @classmethod
    def setUpClass(cls):
        cls.config_fn = os.path.join(TEST_DATA_PATH, 'configs', 'sffp.conf')
        args, configs = parse_args(shlex.split('list -I --config-path {}'.format(cls.config_fn)))
        cls.connection = handlers.OMEROConnection(args, configs)

    @classmethod
    def tearDownClass(cls):
        del cls.connection

    def test_connect(self):
        '''Test that we can connect'''
        with self.connection:
            connected = self.connection.connected
        self.assertTrue(connected)

    def test_images(self):
        '''Test that we can get images'''
        with self.connection:
            images = [image for image in self.connection.images()]
        self.assertGreater(len(images), 0)

    def test_rois(self):
        '''Test that we can get ROIs'''
        with self.connection:
            rois = [roi for roi in self.connection.rois()]
        self.assertGreater(len(rois), 0)

    def test_projects(self):
        '''Test that we can get available projects'''
        with self.connection:
            projects = [project for project in self.connection.projects]
        self.assertGreater(len(projects), 0)

    def test_datasets(self):
        '''test that we can get available datasets'''
        with self.connection:
            datasets = [dataset for dataset in self.connection.datasets()]
        self.assertGreater(len(datasets), 0)

    def test_get_image_by_id(self):
        '''Test that we can get a single image by ID'''
        with self.connection:
            image = self.connection.getImage(101)
        self.assertIsInstance(image, _ImageWrapper)

    def test_get_rois_by_image_id(self):
        '''Test that we can get all ROIs associated with an image'''
        with self.connection:
            rois = self.connection.getROIs(101)
        self.assertGreater(len(rois), 0)


class TestHandlers_OMEROROI(unittest.TestCase):
    def test_create(self):
        self.assertTrue(False)

    def test_addShape(self):
        self.assertTrue(False)

    def test_load_data(self):
        self.assertTrue(False)


class TestPrimitives_Shape(unittest.TestCase):
    '''
    Test generic Shape attributes with a Rectangle ROI object
        FillColor
        FontFamily
        FontSize
        FontStyle
        StrokeColor
        StrokeWidth
        TheC
        TheT
        TheZ
        TextValue
    '''
    @classmethod
    def setUpClass(cls):
        # the connectiona args
        cls.config_fn = os.path.join(TEST_DATA_PATH, 'configs', 'sffp.conf')
        args, configs = parse_args(shlex.split('list -I --config-path {}'.format(cls.config_fn)))
        # the connection
        cls.connection = handlers.OMEROConnection(args, configs)
        # image
        with cls.connection:
            cls.image = cls.connection.getImage(101)
        # the ROI
        cls.roi = handlers.OMEROROI('my_roi', cls.image)
        # vars
        cls.theT = _random_integer()
        cls.theZ = _random_integer()
        cls.X = _random_integer()
        cls.Y = _random_integer()
        cls.Width = _random_integer()
        cls.Height = _random_integer()
        cls.TextValue = 'I am a rectangle'
        cls.r, cls.g, cls.b = _random_float(), _random_float(), _random_float()
        cls.fontSize = _random_integer(5, 15)

        # create the shape with the above params
        cls.rect = primitives.Rectangle(theT=cls.theT, theZ=cls.theZ, X=cls.X, Y=cls.Y, Width=cls.Width, Height=cls.Height)
        cls.rect.setFillColor(cls.r, cls.g, cls.b)
        cls.rect.setFontFamily('serif')
        cls.rect.setFontSize(cls.fontSize, 'POINT')
        cls.rect.setFontStyle('italic')
        cls.rect.setStrokeColor(cls.r, cls.g, cls.b)
        cls.rect.setTextValue(cls.TextValue)

        # add the shape to the ROI
        cls.roi.addShape(cls.rect)

    def test_added_shape(self):
        self.assertTrue(isinstance(self.roi.getShape(0), RectangleI))

    def test_rectangle_type(self):
        self.assertTrue(isinstance(self.rect, RectangleI))

    def test_theT(self):
        self.assertEqual(self.rect.getTheT(), self.theT)

    def test_theZ(self):
        self.assertEqual(self.rect.getTheZ(), self.theZ)

    def test_X(self):
        self.assertEqual(self.rect.getX(), self.X)

    def test_Y(self):
        self.assertEqual(self.rect.getY(), self.Y)

    def test_Width(self):
        self.assertEqual(self.rect.getWidth(), self.Width)

    def test_Height(self):
        self.assertEqual(self.rect.getHeight(), self.Height)

    def test_fillColor(self):
        self.assertEqual(self.rect.getFillColor(), colorInt(self.r, self.g, self.b, 1))

    def test_fontFamily(self):
        self.assertEqual(self.rect.getFontFamily(), 'serif')

    def test_fontSize(self):
        self.assertEqual(self.rect.getFontSize(), self.fontSize)

    def test_fontStyle(self):
        self.assertEqual(self.rect.getFontStyle(), 'italic')

    def test_strokeColor(self):
        self.assertEqual(self.rect.getStrokeColor(), colorInt(self.r, self.g, self.b, 1))

    def test_textValue(self):
        self.assertEqual(self.rect.getTextValue(), self.TextValue)


class TestPrimitives_Line(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X1 = _random_float()
        cls.X2 = _random_float()
        cls.Y1 = _random_float()
        cls.Y2 = _random_float()
        cls.line = primitives.Line(X1=cls.X1, X2=cls.X2, Y1=cls.Y1, Y2=cls.Y2)

    def test_create(self):
        '''Test that we can create a Line object'''
        # assertions
        self.assertEqual(round(self.line.getX1(), 6), round(self.X1, 6))
        self.assertEqual(round(self.line.getX2(), 6), round(self.X2, 6))
        self.assertEqual(round(self.line.getY1(), 6), round(self.Y1, 6))
        self.assertEqual(round(self.line.getY2(), 6), round(self.Y2, 6))

    def test_modify(self):
        '''Test that we can modify params'''
        X1 = _random_float()
        X2 = _random_float()
        Y1 = _random_float()
        Y2 = _random_float()
        self.line.setX1(X1)
        self.line.setX2(X2)
        self.line.setY1(Y1)
        self.line.setY2(Y2)
        # assertions
        self.assertEqual(round(self.line.getX1(), 6), round(X1, 6))
        self.assertEqual(round(self.line.getX2(), 6), round(X2, 6))
        self.assertEqual(round(self.line.getY1(), 6), round(Y1, 6))
        self.assertEqual(round(self.line.getY2(), 6), round(Y2, 6))


class TestPrimitives_Rectangle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X = _random_float() * _random_integer()
        cls.Y = _random_float() * _random_integer()
        cls.Width = _random_float() * _random_integer()
        cls.Height = _random_float() * _random_integer()
        cls.rectangle = primitives.Rectangle(X=cls.X, Y=cls.Y, Width=cls.Width, Height=cls.Height)

    def test_create(self):
        '''Test that we can create a Rectangle object'''
        self.assertEqual(round(self.rectangle.getX(), 6), round(self.X, 6))
        self.assertEqual(round(self.rectangle.getY(), 6), round(self.Y, 6))
        self.assertEqual(round(self.rectangle.getWidth(), 6), round(self.Width, 6))
        self.assertEqual(round(self.rectangle.getHeight(), 6), round(self.Height, 6))

    def test_modify(self):
        '''Test that we can modify params'''
        X = _random_float() * _random_integer()
        Y = _random_float() * _random_integer()
        Width = _random_float() * _random_integer()
        Height = _random_float() * _random_integer()
        self.rectangle.setX(X)
        self.rectangle.setY(Y)
        self.rectangle.setWidth(Width)
        self.rectangle.setHeight(Height)
        # assertions
        self.assertEqual(round(self.rectangle.getX(), 6), round(X, 6))
        self.assertEqual(round(self.rectangle.getY(), 6), round(Y, 6))
        self.assertEqual(round(self.rectangle.getWidth(), 6), round(Width, 6))
        self.assertEqual(round(self.rectangle.getHeight(), 6), round(Height, 6))


class TestPrimitives_Mask(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X = _random_float() * _random_integer()
        cls.Y = _random_float() * _random_integer()
        cls.Width = _random_float() * _random_integer()
        cls.Height = _random_float() * _random_integer()
        cls.mask = primitives.Mask(X=cls.X, Y=cls.Y, Width=cls.Width, Height=cls.Height)

    def test_create(self):
        '''Test that we can create a Mask object'''
        self.assertEqual(round(self.mask.getX(), 6), round(self.X, 6))
        self.assertEqual(round(self.mask.getY(), 6), round(self.Y, 6))
        self.assertEqual(round(self.mask.getWidth(), 6), round(self.Width, 6))
        self.assertEqual(round(self.mask.getHeight(), 6), round(self.Height, 6))

    def test_modify(self):
        '''Test that we can modify params'''
        X = _random_float() * _random_integer()
        Y = _random_float() * _random_integer()
        Width = _random_float() * _random_integer()
        Height = _random_float() * _random_integer()
        self.mask.setX(X)
        self.mask.setY(Y)
        self.mask.setWidth(Width)
        self.mask.setHeight(Height)
        # assertions
        self.assertEqual(round(self.mask.getX(), 6), round(X, 6))
        self.assertEqual(round(self.mask.getY(), 6), round(Y, 6))
        self.assertEqual(round(self.mask.getWidth(), 6), round(Width, 6))
        self.assertEqual(round(self.mask.getHeight(), 6), round(Height, 6))


class TestPrimitives_Ellipse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X = _random_float() * _random_integer()
        cls.Y = _random_float() * _random_integer()
        cls.RadiusX = _random_float() * _random_integer()
        cls.RadiusY = _random_float() * _random_integer()
        cls.ellipse = primitives.Ellipse(X=cls.X, Y=cls.Y, RadiusX=cls.RadiusX, RadiusY=cls.RadiusY)

    def test_create(self):
        '''Test that we can create an Ellipse object'''
        self.assertEqual(round(self.ellipse.getX(), 6), round(self.X, 6))
        self.assertEqual(round(self.ellipse.getY(), 6), round(self.Y, 6))
        self.assertEqual(round(self.ellipse.getRadiusX(), 6), round(self.RadiusX, 6))
        self.assertEqual(round(self.ellipse.getRadiusY(), 6), round(self.RadiusY, 6))

    def test_modify(self):
        '''Test that we can modify params'''
        X = _random_float() * _random_integer()
        Y = _random_float() * _random_integer()
        RadiusX = _random_float() * _random_integer()
        RadiusY = _random_float() * _random_integer()
        self.ellipse.setX(X)
        self.ellipse.setY(Y)
        self.ellipse.setRadiusX(RadiusX)
        self.ellipse.setRadiusY(RadiusY)
        # assertions
        self.assertEqual(round(self.ellipse.getX(), 6), round(X, 6))
        self.assertEqual(round(self.ellipse.getY(), 6), round(Y, 6))
        self.assertEqual(round(self.ellipse.getRadiusX(), 6), round(RadiusX, 6))
        self.assertEqual(round(self.ellipse.getRadiusY(), 6), round(RadiusY, 6))


class TestPrimitives_Point(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X = _random_float() * _random_integer()
        cls.Y = _random_float() * _random_integer()
        cls.point = primitives.Point(X=cls.X, Y=cls.Y)

    def test_create(self):
        '''Test that we can create a Point object'''
        self.assertEqual(round(self.point.getX(), 6), round(self.X, 6))
        self.assertEqual(round(self.point.getY(), 6), round(self.Y, 6))

    def test_modify(self):
        '''Test that we can modify params'''
        X = _random_float() * _random_integer()
        Y = _random_float() * _random_integer()
        self.point.setX(X)
        self.point.setY(Y)
        # assertions
        self.assertEqual(round(self.point.getX(), 6), round(X, 6))
        self.assertEqual(round(self.point.getY(), 6), round(Y, 6))


class TestPrimitives_Polygon(unittest.TestCase):
    '''Tests both primitives.Polyline and primitives.Polygon'''
    @classmethod
    def setUpClass(cls):
        cls.points = [(_random_float() * _random_integer(), _random_float() * _random_integer()) for _ in range(4)]
        cls.polygon = primitives.Polygon()

    def test_create(self):
        '''Test that we can create a Polygon object'''
        self.assertEqual(len(self.polygon.getPoints()), 0)

    def test_modify(self):
        '''Test that we can modify params'''
        self.polygon.setPoints(self.points)
        # assertions
        self.assertEqual(len(self.polygon.getPoints()), 4)

    def test_swapXY(self):
        '''Test that swapping X and Y works'''
        self.polygon.setPoints(self.points, swapXY=True)
        swapped_points = map(lambda x: (x[1], x[0]), self.points)
        self.assertItemsEqual(map(lambda p: (round(p[0], 6), round(p[1], 6)), self.polygon.getPoints()), map(lambda p: (round(p[0], 6), round(p[1], 6)), swapped_points))

    def test_offsetXFrom(self):
        '''Test that we can define an X offset for point values'''
        X_offset = _random_integer()
        self.polygon.setPoints(self.points, offsetXFrom=X_offset)
        offsetX_points = map(lambda p: (X_offset - p[0], p[1]), self.points)
        self.assertItemsEqual(map(lambda p: (round(p[0], 6), round(p[1], 6)), self.polygon.getPoints()), map(lambda p: (round(p[0], 6), round(p[1], 6)), offsetX_points))

    def test_offsetYFrom(self):
        '''Test that we can define an Y offset for point values'''
        Y_offset = _random_integer()
        self.polygon.setPoints(self.points, offsetYFrom=Y_offset)
        offsetY_points = map(lambda p: (p[0], Y_offset - p[1]), self.points)
        self.assertItemsEqual(map(lambda p: (round(p[0], 6), round(p[1], 6)), self.polygon.getPoints()), map(lambda p: (round(p[0], 6), round(p[1], 6)), offsetY_points))


class TestPrimitives_Label(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X = _random_float() * _random_integer()
        cls.Y = _random_float() * _random_integer()
        cls.label_text = "Some very interesting label"
        cls.label = primitives.Label(X=cls.X, Y=cls.Y, label=cls.label_text)

    def test_create(self):
        '''Test that we can create a Point object'''
        self.assertEqual(round(self.label.getX(), 6), round(self.X, 6))
        self.assertEqual(round(self.label.getY(), 6), round(self.Y, 6))
        self.assertEqual(self.label.getTextValue(), self.label_text)

    def test_modify(self):
        '''Test that we can modify params'''
        X = _random_float() * _random_integer()
        Y = _random_float() * _random_integer()
        label_text = "Another piece of label text"
        self.label.setX(X)
        self.label.setY(Y)
        self.label.setTextValue(label_text)
        # assertions
        self.assertEqual(round(self.label.getX(), 6), round(X, 6))
        self.assertEqual(round(self.label.getY(), 6), round(Y, 6))
        self.assertEqual(self.label.getTextValue(), label_text)


class TestHandler_OMEROROI_attachROI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config_fn = os.path.join(TEST_DATA_PATH, 'configs', 'sffp.conf')
        args, configs = parse_args(shlex.split('list -I --config-path {}'.format(cls.config_fn)))
        cls.connection = handlers.OMEROConnection(args, configs)

    @classmethod
    def tearDownClass(cls):
        del cls.connection

    def test_attachRois(self):
        '''Test that we can attach ROIs
          
        Implicitly tests .saveRoi() method
        '''
        from ..formats.roi import ROISegmentation
        args, configs = parse_args(shlex.split('attachroi --config-path {} file.roi'.format(self.config_fn)))
        roi_fn = os.path.join(TEST_DATA_PATH, 'roi', 'test_emd_1832.roi')
        roi_seg = ROISegmentation(roi_fn)
        image_ids = get_image_ids(roi_seg, args)
        # delete rois first
        map(self.delete_rois, image_ids.values())
        rois_before_attach = 0
        rois_after_attach = 0
        with self.connection:
            for roi in self.connection.rois(project='test', dataset='test_attach'):
                rois_before_attach += len(roi[1])
            self.assertEqual(rois_before_attach, 0)
            for orientation in roi_seg.oriented_segments:
                image_id = image_ids[orientation]
                image = self.connection.getImage(image_id)
                omero_rois = roi_seg.as_omero_rois(orientation, image, args)
                self.connection.attachRois(omero_rois)
            for roi in self.connection.rois(project='test', dataset='test_attach'):
                rois_after_attach += len(roi[1])
            self.assertGreater(rois_after_attach, 0)
        # assertions
        # delete rois last
        map(self.delete_rois, image_ids.values())

    def delete_rois(self, image_id):
        '''Delete all ROIs'''
        with self.connection:
            rois = self.connection.getROIs(image_id)
            for roi in rois:
                self.connection.deleteRoi(roi.getId().getValue())
