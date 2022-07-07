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
from . import TEST_DATA_PATH
from ..core.parser import parse_args
from ..formats import roi
from ..schema import SFFPSegmentation
import psycopg2


class TestFormats(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.roi_file = os.path.join(TEST_DATA_PATH, 'roi', 'test_emd_1832.roi')
        cls.sff_file = os.path.join(TEST_DATA_PATH, 'sff', 'test_emd_1832.sff')
        cls.config_fn = os.path.join(TEST_DATA_PATH, 'configs', 'sffp.conf')
        cls.roi_segmentation = roi.ROISegmentation(cls.roi_file)
        cls.sff_segmentation = SFFPSegmentation(cls.sff_file)

    def test_roi_read(self):
        """Test that we can read .roi files"""
        self.assertIsInstance(self.roi_segmentation.header, roi.ROIHeader)
        self.assertEqual(len(self.roi_segmentation.segments), 6)

    def test_from_vtk(self):
        """Test that we can create an ROISegmentation from VTK"""
        args, configs = parse_args(
                'roi create {} -f roi -I emd_1832 '.format(os.path.join(TEST_DATA_PATH, 'sff', 'test_emd_1832_v0.8.0.dev1.sff')), use_shlex=True)
        print(args)
        vtk_segmentation = self.sff_segmentation.as_vtk(args, configs)
        print(vtk_segmentation)
        # to get contours we must first slice!
        vtk_segmentation.slice()
        roi_segmentation = vtk_segmentation.as_roi(args, configs)
        print(roi_segmentation.segments)
        self.assertEqual(len(roi_segmentation.segments), len(self.roi_segmentation.segments))
        self.assertIsInstance(roi_segmentation.segments[0], type(self.roi_segmentation.segments[0]))

    # def test_roi_json(self):
    #     """Test that we can write ROIs as JSONs"""
    #     args, configs = parse_args('createroi {} -o file.json'.format(self.sff_file), use_shlex=True)
    #     vtk_segmentation = self.sff_segmentation.

    # def test_roi_get_image_size(self):
    #     """Test that we can get the image size"""
    #     args, configs = parse_args('createroi --reset-ids {} --config-path {}'.format(self.roi_file, self.config_fn), use_shlex=True)
    #     size = self.roi_segmentation.header.get_image_size(args, configs)
    #     self.assertTupleEqual(size, (64, 64, 64))
