# -*- coding: utf-8 -*-
"""
configs.py
===========

SFFTK configs
"""

__author__  = 'Paul K. Korir, PhD'
__email__   = 'pkorir@ebi.ac.uk, paul.korir@gmail.com'
__date__    = '2016-08-23'


import sys
from collections import OrderedDict
import os.path

from sfftk.core.print_tools import print_date
import sfftkplus


def get_configs(fn=os.path.join(sfftkplus.__path__[0], 'sffplus.conf')):
    """Get SFFTK configs
    
    :param str fn: filename (default: ``./sff.conf``)
    :return dict configs: dictionary of configs
    """ 
    configs = OrderedDict()
    
    try:
        assert os.path.exists(fn)
    except AssertionError:
        raise IOError("File '%s' does not exist." % fn)
    
    with open(fn, 'r') as f:
        for row in f:
            if row[0] == '#': # comments
                continue
            if row.strip() == '': # blank lines
                continue
            name, value = row.strip().split('=')
            configs[name] = value
    
    # assertions
    # CONNECT_WITH must be defined
    try:
        assert 'CONNECT_WITH' in configs
    except AssertionError:
        print_date("Please specify CONNECT_WITH in configs (either LOCAL or REMOTE)")
        sys.exit(1)
    # CONNECT_WITH can only have specified values 
    connect_with_values = ['LOCAL', 'REMOTE']
    try:
        assert configs['CONNECT_WITH'] in connect_with_values
    except AssertionError:
        print_date("Invalid value for var CONNECT_WITH: {}; must be one of: {}".format(configs['CONNECT_WITH'], ", ".join(connect_with_values)))
        sys.exit(1)
    
    return configs


# def set_configs(var, val, fn='./sff.conf'):
#     """Set the variable 'var' to the value held in 'val'
#     
#     :param str var: variable name
#     :param str val: value to set
#     :param str fn: path to the conf file
#     :return int status: 0 on success; otherwise, fail
#     """
#     configs = get_configs(fn=fn)
    