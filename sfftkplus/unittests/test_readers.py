#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sfftkplus.unittests.test_readers

This testing module should have no side-effects because it only reads.
"""

from __future__ import division

import os
import unittest

import __init__ as tests

from ..readers import roireader


__author__  = "Paul K. Korir, PhD"
__email__   = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__    = "2017-08-14"


class TestReaders_roireader(unittest.TestCase):
    def test_get_data(self):
        """Test get_data function"""
        roi_fn = os.path.join(tests.TEST_DATA_PATH, 'test_emd_1832.roi')
        roi_data = roireader.get_data(roi_fn)
        self.assertEqual(roi_data.image_ids.top, 106)
        self.assertEqual(roi_data.image_ids.front, 104)
        self.assertEqual(roi_data.image_ids.right, 105)
        self.assertEqual(len(roi_data.segment), 6)
        self.assertGreater(roi_data.segment[0].colour.red, 0)
        self.assertGreater(roi_data.segment[0].colour.green, 0)
        self.assertGreater(roi_data.segment[0].colour.blue, 0)
        self.assertGreater(roi_data.segment[0].colour.alpha, 0)