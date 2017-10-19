# -*- coding: utf-8 -*-
# schema package
"""
Extended adapter for emdb_sff.py to support ROI handling
"""
__author__  = "Paul K. Korir, PhD"
__email__   = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__    = "2016-09-14"


import sfftk.schema
from sfftkplus.formats.vtkmesh import VTKSegmentation


class SFFPSegmentation(sfftk.schema.SFFSegmentation):
    """Extension of sfftk.schema.SFFSegmentation incorporation ROI conversion"""
    def __init__(self, *args, **kwargs):
        super(SFFPSegmentation, self).__init__(*args, **kwargs)
    def as_vtk(self, args):
        return VTKSegmentation(self, args)