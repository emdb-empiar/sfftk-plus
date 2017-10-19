# -*- coding: utf-8 -*-
# roireader



__author__  = "Paul K. Korir, PhD"
__email__   = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__    = "2017-06-29"


from ..schema.roi import parse


def get_data(fn, *args, **kwargs):
    roi_data = parse(fn, silence=True)
    return roi_data