# -*- coding: utf-8 -*-
"""
configs.py
===========

SFFTK configs
"""

__author__ = 'Paul K. Korir, PhD'
__email__ = 'pkorir@ebi.ac.uk, paul.korir@gmail.com'
__date__ = '2016-08-23'


import sys
import os.path

from sfftk.core.configs import Configs
from sfftk.core.print_tools import print_date
from sfftkplus import BASE_DIR


class SFFPConfigs(Configs):
    shipped_configs = os.path.join(BASE_DIR, 'sffp.conf')

    def read(self):
        super(SFFPConfigs, self).read()
        # CONNECT_WITH must be defined in sffp.conf
        connect_with_values = ['LOCAL', 'REMOTE']
        try:
            assert 'CONNECT_WITH' in self
            # CONNECT_WITH can only have specified values
            try:
                assert self['CONNECT_WITH'] in connect_with_values
            except AssertionError:
                print_date("Invalid value for var CONNECT_WITH: {}; must be one of: {}".format(self['CONNECT_WITH'], ", ".join(connect_with_values)))
        except AssertionError:
            print_date("CONNECT_WITH not found in current configs (with value: {}".format(", ".join(connect_with_values)))
