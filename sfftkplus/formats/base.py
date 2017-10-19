#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
"""
sfftkplus.formats.base



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

__author__  = "Paul K. Korir, PhD"
__email__   = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__    = "2017-04-11"

import sfftk.formats.base


class Mesh(sfftk.formats.base.Mesh):
    pass


class Contours(sfftk.formats.base.Contours):
    pass


class Segment(sfftk.formats.base.Segment):
    pass


class Annotation(sfftk.formats.base.Annotation):
    pass


class Header(sfftk.formats.base.Header):
    pass


class Segmentation(sfftk.formats.base.Segmentation):
    """Extended Segmentation base class"""
    def export(self, *args, **kwargs):
        """Overridden export to handle for ROI"""
        print "Exporting segmentation as a .roi file..."
        