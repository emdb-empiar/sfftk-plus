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
import unittest

from . import TEST_DATA_PATH
from ..core.parser import parse_args
from ..formats import roi
from ..schema import SFFPSegmentation


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
            'roi create {} -f roi -I emd_1832 '.format(
                os.path.join(TEST_DATA_PATH, 'sff', 'test_emd_1832_v0.8.0.dev1.sff')), use_shlex=True)
        print(args)
        vtk_segmentation = self.sff_segmentation.as_vtk(args, configs)
        print(vtk_segmentation)
        # to get contours we must first slice!
        vtk_segmentation.slice()
        roi_segmentation = vtk_segmentation.as_roi(args, configs)
        print(roi_segmentation.segments)
        self.assertEqual(len(roi_segmentation.segments), len(self.roi_segmentation.segments))
        self.assertIsInstance(roi_segmentation.segments[0], type(self.roi_segmentation.segments[0]))

    def test_export(self):
        """Test that we can export as vtp files"""
        args, configs = parse_args(
            'export {}'.format(
                os.path.join(TEST_DATA_PATH, 'sff', 'test_emd_1832_v0.8.0.dev1.json')
            ), use_shlex=True
        )
        sff_seg = SFFPSegmentation.from_file(args.sff_file)
        vtk_segmentation = sff_seg.as_vtk(args, configs)
        vtk_segmentation.export('original_transform', args, configs)

