# -*- coding: utf-8 -*-
# test_schema.py

import os
import shlex
import unittest

import h5py

import __init__ as tests
from ..core.parser import parse_args
from ..formats import vtkmesh
from ..schema import SFFPSegmentation


__author__  = "Paul K. Korir, PhD"
__email__   = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__    = "2017-08-17"


class TestSFFPSegmentation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sff_file = os.path.join(tests.TEST_DATA_PATH, 'sff', 'test_emd_1832.sff')
        cls.hff_file = os.path.join(tests.TEST_DATA_PATH, 'sff',  'test_emd_1832.hff')
        cls.json_file = os.path.join(tests.TEST_DATA_PATH,  'sff', 'test_emd_1832.json')
        
    def test_read_sff(self):
        """Test that we can read an .sff file"""
        sff_segmentation = SFFPSegmentation(self.sff_file)
        self.assertEqual(sff_segmentation.version, '0.7.0.dev0')
        self.assertEqual(len(sff_segmentation.segments), 6)
        
    def test_read_hff(self):
        """Test that we can read an .hff file"""
        with h5py.File(self.hff_file) as h:
            hff_segmentation = SFFPSegmentation.from_hff(h)
        self.assertEqual(hff_segmentation.version, '0.7.0.dev0')
        self.assertEqual(len(hff_segmentation.segments), 6)
        
    def test_read_json(self):
        """Test that we can read a .json file"""
        json_segmentation = SFFPSegmentation.from_json(self.json_file)
        self.assertEqual(json_segmentation.version, '0.7.0.dev0')
        self.assertEqual(len(json_segmentation.segments), 6)
        
    def test_as_vtk(self):
        """Test that we can convert to a VTK object"""
        args, configs = parse_args(shlex.split('createroi file.sff -o file.roi'))
        sff_segmentation = SFFPSegmentation(self.sff_file)
        vtk_segmentation = sff_segmentation.as_vtk(args, configs)
        self.assertIsInstance(vtk_segmentation, vtkmesh.VTKSegmentation)

        