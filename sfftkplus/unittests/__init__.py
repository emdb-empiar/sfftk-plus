# -*- coding: utf-8 -*-

# path to test data
import os
import random

from .. import BASE_DIR


TEST_DATA_PATH = os.path.join(BASE_DIR, 'test_data')

# helper functions
def _random_integer(start=1, stop=1000): return random.randint(start, stop)
def _random_float(): return random.random()


