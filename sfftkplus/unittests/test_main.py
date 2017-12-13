#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

"""
test_py

Unit tests for convert subcommand
"""

__author__  = 'Paul K. Korir, PhD'
__email__   = 'pkorir@ebi.ac.uk'
__date__    = '2016-06-10'


import sys
import os
import glob
import unittest
from ..core.parser import parse_args
import __init__ as tests # import _random_integer, _random_float

sys.path.insert(0, "..")


# redirect sys.stderr/sys.stdout to /dev/null
# from: http://stackoverflow.com/questions/8522689/how-to-temporary-hide-stdout-or-stderr-while-running-a-unittest-in-python
_stderr = sys.stderr
_stdout = sys.stdout
null = open(os.devnull, 'wb')
sys.stdout = null
sys.stderr = null

user = 'test_user'
password = 'test'
host = 'localhost'
port = '4064'


#===============================================================================
# sfftk-plus tests
#===============================================================================
# class TestMain_handle_list(unittest.TestCase):
#     def test_list_images(self):
#         """Test that we can list all images"""
#         cmd = 'list -U {} -P {} -H {} -p {} -I'.format(user, password, host, port ).split(' ')
#         args = parse_args(cmd)
#           
#         self.assertEqual(0, Main.handle_list(args))
#      
#     def test_list_image(self):
#         """Test that we can list one image"""
#         cmd = 'list -U {} -P {} -H {} -p {} -I -i 401'.format(user, password, host, port ).split(' ')
#         args = parse_args(cmd)
#           
#         self.assertEqual(0, Main.handle_list(args))
#      
#     def test_list_rois(self):
#         """Test listing of all ROIs"""
#         cmd = 'list -U {} -P {} -H {} -p {} -R'.format(user, password, host, port ).split(' ')
#         args = parse_args(cmd)
#           
#         self.assertEqual(0, Main.handle_list(args))
#  
#  
# class TestMain_handle_createroi(unittest.TestCase):
#     @staticmethod
#     def clear_files():
#         map(os.remove, glob.glob('sff/test_data/*.sff'))
#         map(os.remove, glob.glob('sff/test_data/*.roi'))
#         
#     def setUp(self):
#         """Clear .sff and .roi files"""
#         unittest.TestCase.setUp(self)
#         self.clear_files()       
#         
#     def tearDown(self):
#         """Clear .sff and .roi files"""
#         unittest.TestCase.tearDown(self)
#         self.clear_files()
#         
#     def test_rois_from_meshlist(self):
#         """Test that we can make an ROI file from an SFF file"""
#         args = parse_args(['convert', '-o', 'sff/test_data/test_data.sff', 'sff/test_data/test_data.surf'])
#          
#         Main.handle_convert(args)
#          
#         sff_files = glob.glob('sff/test_data/*.sff')        
#          
#         self.assertEqual(len(sff_files), 1)
#          
#         args = parse_args('createroi -o sff/test_data/test_data.roi -P meshList sff/test_data/test_data.sff'.split(' '))
#          
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#      
#     def test_rois_from_amira_surf(self):
#         """Test that we can make an ROI file from an Amira .surf file"""
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.surf'.split(' '))
#          
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#      
#     def test_rois_from_imod(self):
#         """Test that we can make an ROI file from an IMOD file"""
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.mod'.split(' '))
#          
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#          
#     def test_rois_from_segger(self):
#         """Test that we can make an ROI file from a Segger file"""
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.seg'.split(' '))
#          
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#     
#     def test_rois_from_emdb_map(self):
#         """Test that we can make an ROI file from an EMDB .map file"""
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.map'.split(' '))
#         
#         Main.handle_createroi(args)
#         
#         roi_files = glob.glob('sff/test_data/*.roi')
#         
#         self.assertEqual(len(roi_files), 1)
#     
#     def test_rois_from_stl(self):
#         """Test that we can make an ROI file from a StL file"""
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.stl'.split(' '))
#          
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#         
#     
# #     def test_rois_fail_from_contourlist(self):
# #         """Test that we raise an exception with a contourList"""
# #         args = parse_args('convert -o sff/test_data/test_data.sff -d contourList sff/test_data/test_data.mod'.split(' '))
# #          
# #         Main.handle_convert(args)
# #          
# #         sff_files = glob.glob('sff/test_data/*.sff')        
# #          
# #         self.assertEqual(len(sff_files), 1)
# #          
# #         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.sff'.split(' '))
# #          
# #         with self.assertRaises(ValueError):
# #             Main.handle_createroi(args)
# 
# 
# class TestMain_handle_attachroi(unittest.TestCase):
#     @staticmethod
#     def clear_files():
#         map(os.remove, glob.glob('sff/test_data/*.sff'))
#         map(os.remove, glob.glob('sff/test_data/*.roi'))
#          
#     def setUp(self):
#         """Clear .sff and .roi files"""
#         unittest.TestCase.setUp(self)
#         self.clear_files()       
#          
#     def tearDown(self):
#         """Clear .sff and .roi files"""
#         unittest.TestCase.tearDown(self)
#         self.clear_files()
#          
#     def test_attach_success(self):
#         """Test successful attach of ROIs to image"""
#         args = parse_args('convert -o sff/test_data/test_data.sff sff/test_data/test_data.mod'.split(' '))
#          
#         Main.handle_convert(args)
#          
#         sff_files = glob.glob('sff/test_data/*.sff')
#                  
#         self.assertEqual(len(sff_files), 1)
#          
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.sff'.split(' '))
#         
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#          
#         args = parse_args('attachroi -U {} -P {} -H {} -p {} -x 501 -y 502 -z 503 sff/test_data/test_data.roi'.format(user, password, host, port ).split(' '))
#          
#         self.assertEqual(0, Main.handle_attachroi(args))
#     
#     def test_attach_fail(self):
#         """Test that attach will fail for the wrong image ID"""
#         args = parse_args('convert -o sff/test_data/test_data.sff sff/test_data/test_data.mod'.split(' '))
#          
#         Main.handle_convert(args)
#          
#         sff_files = glob.glob('sff/test_data/*.sff')
#          
#         self.assertEqual(len(sff_files), 1)
#          
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.sff'.split(' '))
#         
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#         
#         args = parse_args('attachroi -U {} -P {} -H {} -p {} -x 10 sff/test_data/test_data.roi'.format(user, password, host, port).split(' '))
#         
#         self.assertEqual(1, Main.handle_attachroi(args))  
#             
# 
# class TestMain_handle_deleteroi(unittest.TestCase):
#     @staticmethod
#     def clear_files():
#         map(os.remove, glob.glob('sff/test_data/*.sff'))
#         map(os.remove, glob.glob('sff/test_data/*.roi'))
#          
#     def setUp(self):
#         """Clear .sff and .roi files"""
#         unittest.TestCase.setUp(self)
#         self.clear_files()       
#          
#     def tearDown(self):
#         """Clear .sff and .roi files"""
#         unittest.TestCase.tearDown(self)
#         self.clear_files()
#         
#     def test_delete_single(self):
#         """Test successful delete of one ROI on image"""
#         args = parse_args('convert -o sff/test_data/test_data.sff sff/test_data/test_data.mod'.split(' '))
#          
#         Main.handle_convert(args)
#          
#         sff_files = glob.glob('sff/test_data/*.sff')
#          
#         self.assertEqual(len(sff_files), 1)
#          
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.sff'.split(' '))
#         
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#          
#         args = parse_args('attachroi -U {} -P {} -H {} -p {} -x 501 -y 502 -z 503 sff/test_data/test_data.roi'.format(user, password, host, port ).split(' '))
#          
#         self.assertEqual(0, Main.handle_attachroi(args))
#         
#         from omero_wrapper.handlers import OMEROConnectionHandler
#         import random
#         
#         # get an roi ID
#         with OMEROConnectionHandler(user=args.user, password=args.password, host=args.host, port=args.port) as omero_connection:
#             roi = random.choice(omero_connection.rois(501))
#             roi_id = roi.getId().getValue()
#             
#             
#             args = parse_args('deleteroi -U {} -P {} -H {} -p {} -r {}'.format(user, password, host, port , roi_id).split(' '))
#         
#             self.assertEqual(0, Main.handle_deleteroi(args))
#             
#             roi_ids = map(lambda x: x.getId().getValue(), omero_connection.rois(501))
#             
#             self.assertNotIn(roi_id, roi_ids)
#     
#     def test_delete_multiple(self):
#         """Test successful delete of all ROIs on image"""
#         args = parse_args('convert -o sff/test_data/test_data.sff sff/test_data/test_data.mod'.split(' '))
#          
#         Main.handle_convert(args)
#          
#         sff_files = glob.glob('sff/test_data/*.sff')
#          
#         self.assertEqual(len(sff_files), 1)
#          
#         args = parse_args('createroi -o sff/test_data/test_data.roi sff/test_data/test_data.sff'.split(' '))
#         
#         Main.handle_createroi(args)
#          
#         roi_files = glob.glob('sff/test_data/*.roi')
#          
#         self.assertEqual(len(roi_files), 1)
#          
#         args = parse_args('attachroi -U {} -P {} -H {} -p {} -x 501 -y 502 -z 503 sff/test_data/test_data.roi'.format(user, password, host, port ).split(' '))
#          
#         self.assertEqual(0, Main.handle_attachroi(args))
#         
#         from omero_wrapper.handlers import OMEROConnectionHandler
#         import random
#         
#         # get an roi ID
#         with OMEROConnectionHandler(user=args.user, password=args.password, host=args.host, port=args.port) as omero_connection:
#             args = parse_args('deleteroi -U {} -P {} -H {} -p {} -i {}'.format(user, password, host, port , 501).split(' '))
#         
#             self.assertEqual(0, Main.handle_deleteroi(args))
#             
#             roi_ids = map(lambda x: x.getId().getValue(), omero_connection.rois(501))
#             
#             self.assertEqual(0, len(roi_ids))
#             
#     
# class TestMain_handle_view3d(unittest.TestCase):
#     def test_renderer(self):
#         """Test that we can render a mesh from an .sff file
#         
#         The test consists of seeing that the mesh is not empty (it has points and cells)
#         """
#         # read an sff file
#         from readers import sffreader
#         
#         descriptor, segments, colours, alphas = sffreader.get_data('sff/sff_data/test_meshList.sff')
#         # create a mesh
#         from converters import meshtools
#         
#         mesh = meshtools.get_meshes(segments, colours, alphas, form=descriptor)[0]
#         
#         # check that the mesh has points and cells
#         self.assertTrue(mesh.GetNumberOfPoints() > 0)
#         self.assertTrue(mesh.GetNumberOfCells() > 0)
# 
# 
# class TestMain_handle_notes_readonly(unittest.TestCase):   
#     #===========================================================================
#     # notes: list
#     #===========================================================================
#     def test_notes_list(self):
#         """Test that we can list notes"""
#         args = parse_args('notes list sff/sff_data/test_annotations.sff'.split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     def test_notes_list_short(self):
#         """Test that we can short list notes"""
#         args = parse_args('notes list -s sff/sff_data/test_annotations.sff'.split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     #===========================================================================
#     # notes: show
#     #===========================================================================
#     def test_notes_show(self):
#         """Test that we can show notes"""
#         args =  parse_args('notes show -i 104 sff/sff_data/test_annotations.sff'.split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     def test_notes_show_sort(self):
#         """Test that we can show short notes for a subset of IDs"""
#         args = parse_args('notes show -i 104,105 sff/sff_data/test_annotations.sff'.split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     #===========================================================================
#     # notes: search
#     #===========================================================================
#     def test_notes_search_default(self):
#         """Default search without options"""
#         args = parse_args("notes search 'mitochondria'".split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     def test_notes_search_exact(self):
#         """Test exact searches"""
#         args = parse_args("notes search -x mitochondria".split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     def test_notes_search_start(self):
#         """Test that we can search starting at a particular result"""
#         start = _random_integer()
#         args = parse_args("notes search -s {} mitochondria".format(start).split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     def test_notes_search_rows(self):
#         """Test that we can return a certain number of rows"""
#         rows = _random_integer()
#         args = parse_args("notes search -r {} mitochondria".format(rows).split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     def test_notes_search_ontology(self):
#         """Test that we can search a particular ontology"""
#         args = parse_args("notes search -O fma mitochondria".split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
#     def test_notes_search_ontology_fails(self):
#         """Test that we can catch a wrong ontology"""
#         args = parse_args("notes search -O sldkjfl mitochondria".split(' '))
#         self.assertEqual(1, Main.handle_notes(args))
#     def test_notes_search_obsoletes(self):
#         """Test that we can search for obsoletes"""
#         args = parse_args("notes search -o mitochondria".split(' '))
#         self.assertEqual(0, Main.handle_notes(args))
# 
# 
# class TestMain_handle_notes_rw(unittest.TestCase):
#     def setUp(self):
#         unittest.TestCase.setUp(self)
#     def tearDown(self):
#         unittest.TestCase.tearDown(self)
#     def test_notes_add(self):
#         pass
#     def test_notes_edit(self):
#         pass
#     def test_notes_del(self):
#         pass
#     def test_notes_save(self):
#         pass
#     def test_notes_trash(self):
#         pass
    
