# -*- coding: utf-8 -*-
# test_formats.py
"""
sfftkplus.formats modules unit tests
"""
from __future__ import division, print_function

__author__ = "Paul K. Korir, PhD"
__email__ = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__ = "2017-08-14"

import os
import sys
import shlex
import unittest
import __init__ as tests
from ..core.parser import parse_args
from ..formats import roi
from ..schema import SFFPSegmentation
import psycopg2


class TestFormats(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.roi_file = os.path.join(tests.TEST_DATA_PATH, 'roi', 'test_emd_1832.roi')
        cls.sff_file = os.path.join(tests.TEST_DATA_PATH, 'sff', 'test_emd_1832.sff')
        cls.config_fn = os.path.join(tests.TEST_DATA_PATH, 'configs', 'sffp.conf')
        cls.roi_segmentation = roi.ROISegmentation(cls.roi_file)
        cls.sff_segmentation = SFFPSegmentation(cls.sff_file)

    def test_roi_read(self):
        """Test that we can read .roi files"""
        self.assertIsInstance(self.roi_segmentation.header, roi.ROIHeader)
        self.assertEqual(len(self.roi_segmentation.segments), 6)

    def test_from_vtk(self):
        """Test that we can create an ROISegmentation from VTK"""
        args, configs = parse_args(
            shlex.split(
                'createroi {} -o file.roi'.format(os.path.join(tests.TEST_DATA_PATH, 'sff', 'test_emd_1832.sff'))))
        vtk_segmentation = self.sff_segmentation.as_vtk(args, configs)
        #  to get contours we must first slice!
        vtk_segmentation.slice()
        roi_segmentation = vtk_segmentation.as_roi(args, configs)
        self.assertEqual(len(roi_segmentation.segments), len(self.roi_segmentation.segments))
        self.assertIsInstance(roi_segmentation.segments[0], type(self.roi_segmentation.segments[0]))

    # def test_roi_json(self):
    #     """Test that we can write ROIs as JSONs"""
    #     args, configs = parse_args(shlex.split('createroi {} -o file.json'.format(self.sff_file)))
    #     # vtk_segmentation = self.sff_segmentation.

    def test_roi_get_image_size(self):
        """Test that we can get the image size"""
        args, configs = parse_args(shlex.split('createroi --reset-ids {} --config-path {}'.format(self.roi_file, self.config_fn)))
        size = self.roi_segmentation.header.get_image_size(args, configs)
        self.assertTupleEqual(size, (64, 64, 64))