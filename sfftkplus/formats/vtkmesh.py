#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sfftkplus.formats.vtkmesh



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

from __future__ import division

import math
import os
from random import random

import h5py
import numpy
from sfftk.core.print_tools import print_date
from sfftk.readers.segreader import get_root
from sfftkplus.formats.base import Segmentation, Header, Segment, Contours, Mesh
import vtk


__author__ = "Paul K. Korir, PhD"
__email__ = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__ = "2017-04-12"
__updated__ = "2018-02-27"


def simplify_mask(mask, r_ids, r_p_zip, replace=True):
    """Simplify the mask by replacing all `region_ids` with their `root_parent_id`
    
    The `region_ids` and `parent_ids` are paired from which a tree is inferred. The root 
    of this tree is value `0`. `region_ids` that have a corresponding `parent_id` of 0
    are penultimate roots. This method replaces each `region_id` with its penultimate `parent_id`.
    It *simplifies* the volume.
    
    :param mask: a 3D volume
    :type mask: `numpy.array`
    :param r_id: sequence of `region_id`
    :type r_id: iterable
    :param r_p_zip: sequence of 2-tuples with `region_id` and `parent_id`
    :type r_p_zip: iterable
    :param bool replace: if `True` then the returned `mask` will have values; `False` will leave the `mask` unchanged (useful for running tests to speed things up)
    :return: `simplified_mask`, `segment_colours`, `segment_ids`
    :rtype: tuple
    """
    simplified_mask = numpy.ndarray(mask.shape, dtype=int)  # @UnusedVariable
    simplified_mask = 0

    # group regions_ids by parent_id
    root_parent_id_group = dict()
    for r in r_ids:
        p = get_root(r_p_zip, r)
        if p not in root_parent_id_group:
            root_parent_id_group[p] = [r]
        else:
            root_parent_id_group[p] += [r]

    if replace:
        # It is vastly faster to use multiple array-wide comparisons than to do
        # comparisons element-wise. Therefore, we generate a string to be executed
        #  that will do hundreds of array-wide comparisons at a time.
        # Each comparison is for all region_ids for a parent_id which will
        # then get assigned the parent_id.
        for parent_id, region_id_list in root_parent_id_group.items():
            # check whether any element in the mask has a value == r0 OR r1 ... OR rN
            # e.g. (mask == r0) | (mask == r1) | ... | (mask == rN)
            comp = ' | '.join(['( mask == %s )' % r for r in region_id_list])
            # set those that satisfy the above to have the parent_id
            # Because parent_ids are non-overlapping (i.e. no region_id has two parent_ids)
            # we can do successive summation instead of assignments.
            full_op = 'simplified_mask += (' + comp + ') * %s' % parent_id
            exec(full_op)
    else:
        simplified_mask = mask

    segment_ids = root_parent_id_group.keys()

#     segment_colors = [r_c_zip[s] for s in segment_ids]

    return simplified_mask, segment_ids


class VTKMesh(vtk.vtkPolyData):
    def __init__(self, colour, args):
        if len(colour) == 3:
            self._colour = colour
            self._alpha = 1
        elif len(colour) == 4:
            self._colour, self._alpha = colour[:3], colour[-1]
        self._vtk_args = args

    @property
    def vtk_args(self):
        return self._vtk_args

    @property
    def colour(self):
        return self._colour

    @property
    def alpha(self):
        return self._alpha

    @property
    def edge_colour(self):
        return 0.8, 0.8, 1

    @property
    def smooth(self):
        return self._vtk_args.smooth

    @property
    def smooth_iterations(self):
        return self._vtk_args.smooth_iterations

    @property
    def transparency(self):
        return self._vtk_args.transparency

    @classmethod
    def from_mesh(cls, sff_mesh, colour, args):
        """Initialiase a VTKMesh object from a ``sfftk.schema.SFFMesh``
        
        :param mesh: a mesh with vertices and polygons
        :type mesh: ``sfftk.schema.SFFMesh``
        :param colour: the segment colour
        :type colour: ``sfftk.schema.SFFColour``
        :param args: parsed arguments
        :type args: ``argparse.Namespace``
        :return vtkmesh: an VTKMesh object
        :rtype vtkmesh: ``VTKMesh``  
        """
        # for each vertex in this mesh
        vertices = dict()
        vertex_ids = list()
        normals = dict()
        normal_ids = list()
        for v in sff_mesh.vertices:
            if v.designation == 'vertex':
                vertices[v.vID] = v.point
                vertex_ids.append(v.vID)
            elif v.designation == 'normal':
                normals[v.vID] = v.point
                normal_ids.append(v.vID)

        if len(normal_ids) == len(vertex_ids):
            normal_to_vertex = dict(zip(normal_ids, vertex_ids))
        else:
            normal_to_vertex = dict()
        # for each polygon in this mesh
        polygons = dict()
        for P in sff_mesh.polygons:
            if len(normal_ids) > 0:
                if P.vertex_ids[0] in normal_ids:  # if the first vertex in the polygon is a normal
                    polygons[P.PID] = tuple(P.vertex_ids[1::2])
                elif P.vertex_ids[0] in vertex_ids:  # if the first vertex in the polygon is a vertex
                    polygons[P.PID] = tuple(P.vertex_ids[::2])
            else:
                polygons[P.PID] = tuple(P.vertex_ids)
        """
        The problem with the above process is that in each mesh we end up with indices offset
        by 2 i.e. e.g. vertices will be 0, 2, 4, ... and normal 1, 3, 5, ...
        because they are defined that way in some files. We need to fix that
         
        We need 
        vertices 0, 1, 2,...
        polygons refer to these vertices
        normals 0, 1, 2, refer to these vertices
        then we do away with normal_to_vertex
        """
        new_vertex_id = 0
        new_vertices = dict()
        old_vertex_id_to_new_vertex_id = dict()
        for vertex_id, vertex in vertices.iteritems():
            new_vertices[new_vertex_id] = vertex
            old_vertex_id_to_new_vertex_id[vertex_id] = new_vertex_id
            new_vertex_id += 1

        new_polygons = dict()
        for polygon_id, polygon in polygons.iteritems():
            new_polygons[polygon_id] = tuple([old_vertex_id_to_new_vertex_id[p] for p in polygon])

        new_normals = dict()
        for normal_id, normal in normals.iteritems():
            old_vertex_id = normal_to_vertex[normal_id]
            new_vertex_id = old_vertex_id_to_new_vertex_id[old_vertex_id]
            new_normals[new_vertex_id] = normal
        # the VTKMesh
        vtkmesh = cls(colour, args)
        # define the geometry
        points = vtk.vtkPoints()
        for vertex_id, vertex in new_vertices.iteritems():
            points.InsertPoint(vertex_id, *vertex)
        vtkmesh.SetPoints(points)
        # define the topology
        cellArray = vtk.vtkCellArray()
        for polygon in new_polygons.itervalues():
            cell_size = len(polygon)
            cellArray.InsertNextCell(cell_size, polygon)
        vtkmesh.SetPolys(cellArray)
        #  define normals (if they exist)
        if not args.normals_off:
            if len(new_normals) == len(new_vertices):
                normals = vtk.vtkFloatArray()
                normals.SetNumberOfComponents(3)
                for normal_id, normal in new_normals.iteritems():
                    normals.InsertTuple3(normal_id, *normal)
                vtkmesh.GetPointData().SetNormals(normals)
        return vtkmesh

    @classmethod
    def from_volume(cls, mask, colour, args):
        """Initialiase a VTKMesh object from a binarised ``numpy.ndarray``
        
        :param mask: a mask (3D matrix)
        :type mask: ``numpy.ndarray``
        :param colour: the segment colour
        :type colour: ``sfftk.schema.SFFColour``
        :param args: parsed arguments
        :type args: ``argparse.Namespace``
        :return vtkmesh: an VTKMesh object
        :rtype vtkmesh: ``VTKMesh``  
        """
        vtkmesh = cls(colour, args)
        # try and figure out the mask value if not specified
        if args.mask_value is not None:
            mask_value = args.mask_value
        else:
            mask_values = filter(lambda x: x != 0, set(mask.flatten().tolist()))
            if len(mask_values) == 1:
                mask_value = mask_values[0]
            else:
                # use the mean
                mask_value = sum(mask_values) / len(mask_values)

        vol = vtk.vtkImageData()
        z_size, y_size, x_size = mask.shape
        vol.SetExtent(0, x_size - 1, 0, y_size - 1 , 0, z_size - 1)
        vtkmesh.vol_extent = 0, x_size - 1, 0, y_size - 1 , 0, z_size - 1
        vol.SetOrigin(0.0, 0.0, 0.0)
        sp_x = 1.0  # /(x_size - 1)
        sp_y = 1.0  # /(y_size - 1)
        sp_z = 1.0  # /(z_size - 1)
        vol.SetSpacing(sp_x, sp_y, sp_z)
        vol.AllocateScalars(vtk.VTK_FLOAT, 1)
        voxels = vtk.vtkFloatArray()
        voxels.SetNumberOfComponents(1)
        for m in mask.flatten().tolist():
            voxels.InsertNextValue(float(m))
        vol.GetPointData().SetScalars(voxels)
        # convert the volume to a surface
        contours = vtk.vtkContourFilter()
        contours.SetInputData(vol)
        contours.SetValue(0, mask_value)
        contours.Update()
        contoursOutput = contours.GetOutput()

        vtkmesh.SetPoints(contoursOutput.GetPoints())
        vtkmesh.SetPolys(contoursOutput.GetPolys())
        # triangulate
        triangleMesh = vtk.vtkTriangleFilter()
        triangleMesh.SetInputData(contoursOutput)
        triangleMesh.Update()
        triangleMeshOutput = triangleMesh.GetOutput()
        # decimate
        decimateMesh = vtk.vtkDecimatePro()
        decimateMesh.SetInputData(triangleMeshOutput)
        decimateMesh.SetTargetReduction(0.9)
        decimateMesh.PreserveTopologyOn()
        decimateMesh.Update()
        decimateMeshOutput = decimateMesh.GetOutput()
        # smooth
        smoothedMesh = vtk.vtkSmoothPolyDataFilter()
        smoothedMesh.SetInputData(decimateMeshOutput)
        smoothedMesh.SetNumberOfIterations(200)
        smoothedMesh.Update()
        smoothedMeshOutput = smoothedMesh.GetOutput()
        # triangle strips
        triangleStrips = vtk.vtkStripper()
        triangleStrips.SetInputData(smoothedMeshOutput)
        triangleStrips.SetMaximumLength(1000)
        triangleStrips.Update()
        triangleStripsOutput = triangleStrips.GetOutput()
        # finally set things
        if not args.normals_off:
            # normals
            normals = vtk.vtkPolyDataNormals()
            normals.SetInputData(triangleStripsOutput)
    #         normals.ComputeCellNormalsOn()
            normals.ConsistencyOn()
            normals.AutoOrientNormalsOn()
            normals.ComputePointNormalsOn()
            normals.SplittingOff()
            normals.SetFeatureAngle(240.0)
            normals.Update()
            normalsOutput = normals.GetOutput()
#             print vtkmesh_n.GetPointData().GetNormals()
            vtkmesh.SetPoints(normals.GetOutput().GetPoints())
            vtkmesh.SetPolys(normals.GetOutput().GetPolys())
            vtkmesh.GetPointData().SetNormals(normalsOutput.GetPointData().GetNormals())
            return vtkmesh
        else:
            vtkmesh.SetPoints(triangleStripsOutput.GetPoints())
            vtkmesh.SetPolys(triangleStripsOutput.GetPolys())
            return vtkmesh

    @classmethod
    def from_shape(cls, shape, colour, args, transform, resolution=20):
        """Initialiase a VTKMesh object from a sfftk.schema.SFFShape
        
        :param shapes: an iterable of shapes
        :type shapes: ``sfftk.schema.SFFShapePrimitiveList
        :param colour: the segment colour
        :type colour: ``sfftk.schema.SFFColour``
        :param args: parsed arguments
        :type args: ``argparse.Namespace``
        :param transform: transform bearing this shape's translation from the origin
        :type transform: ``sfftk.schema.SFFTransform`` subclass (typically ``sfftk.schema.SFFTransformationMatrix``
        :param int resolution: mesh resolution
        :return vtkmesh: an VTKMesh object
        :rtype vtkmesh: ``VTKMesh``  
        """
        assert resolution > 0
        vtkmesh = cls(colour, args)
        from sfftk.schema import SFFEllipsoid, SFFCuboid, SFFCylinder, SFFCone
        if isinstance(shape, SFFEllipsoid):
            vtk_shape = vtk.vtkSphereSource()
            vtk_shape.SetRadius(shape.x)
            """
            :TODO: make this variable
            """
            vtk_shape.SetPhiResolution(resolution)
            vtk_shape.SetThetaResolution(resolution)
        elif isinstance(shape, SFFCylinder):
            vtk_shape = vtk.vtkCylinderSource()
            vtk_shape.SetHeight(shape.height)
            vtk_shape.SetRadius(shape.diameter / 2)
            vtk_shape.SetResolution(resolution)
        elif isinstance(shape, SFFCone):
            vtk_shape = vtk.vtkConeSource()
            vtk_shape.SetHeight(shape.height)
            vtk_shape.SetRadius(shape.bottomRadius)
            vtk_shape.SetResolution(resolution)
        elif isinstance(shape, SFFCuboid):
            vtk_shape = vtk.vtkCubeSource()
            vtk_shape.SetXLength(shape.x)
            vtk_shape.SetYLength(shape.y)
            vtk_shape.SetZLength(shape.z)
        T = transform.data_array
        vtk_shape.SetCenter(float(T[0, 3]), float(T[1, 3]), float(T[2, 3]))
        vtk_shape.Update()
        _vtkmesh = vtk_shape.GetOutput()
        # triangle filter
        triangleMesh = vtk.vtkTriangleFilter()
        triangleMesh.SetInputData(_vtkmesh)
        triangleMesh.Update()
        triangleMeshOutput = triangleMesh.GetOutput()
        vtkmesh.SetPoints(triangleMeshOutput.GetPoints())
        vtkmesh.SetPolys(triangleMeshOutput.GetPolys())
        return vtkmesh

    def slice(self):
        return VTKContours(self)

    def translate(self, to):
        # the transform
        translation = vtk.vtkTransform()
        translation.Translate(*to)
        # transform filter
        transformFilter = vtk.vtkTransformFilter()
        transformFilter.SetInputData(self)
        transformFilter.SetTransform(translation)
        transformFilter.Update()
        return transformFilter.GetOutput()


    def render(self, renderer):
        """
        Render the mesh
         
        :param renderer: a renderer into which an actor derived from this mesh will be added
        :type renderer: `vtk.vtkRenderer`
        """
        assert isinstance(renderer, vtk.vtkRenderer)

        self.mapper = vtk.vtkOpenGLPolyDataMapper()
        # normals
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(self)
        normals.ComputePointNormalsOn()
        normals.SplittingOff()
        normals.AutoOrientNormalsOn()
        normals.ConsistencyOn()
        normals.Update()
        self.mapper.SetInputData(normals.GetOutput())
        self.actor = vtk.vtkOpenGLActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(self.colour)
        if self._vtk_args.wireframe:
            self.actor.GetProperty().SetRepresentationToWireframe()
        if self.alpha:
            self.actor.GetProperty().SetOpacity(self.alpha)
        if self._vtk_args.view_edges:
            self.actor.GetProperty().EdgeVisibilityOn()
            self.actor.GetProperty().SetEdgeColor(*self.edge_colour)
        # lighting and shading
        self.actor.GetProperty().SetInterpolationToPhong()
        self.actor.GetProperty().SetDiffuse(0.7)
        self.actor.GetProperty().SetSpecular(0.5)
        self.actor.GetProperty().SetSpecularPower(40)
#         self.actor.GetProperty().BackfaceCullingOn()
#         self.actor.GetProperty().FrontfaceCullingOn()
        renderer.AddActor(self.actor)
        return renderer


class VTKMeshes(Mesh):
    def __init__(self, sff_segment, colour, args, *args_, **kwargs_):
        self._sff_segment = sff_segment
        self._vtk_args = args
        self._vtk_meshes = list()
        self._colour = colour
        # transforms
        if 'transforms' in kwargs_:
            self._transforms = kwargs_['transforms']
        if 'lattice' in kwargs_:
            self._lattice = kwargs_['lattice']
        if self.primary_descriptor == "meshList":
            for mesh in self._sff_segment.meshes:
                self._vtk_meshes.append(
                    VTKMesh.from_mesh(mesh, self.colour, self._vtk_args)
                    )
        elif self.primary_descriptor == "shapePrimitiveList":
            for shape in self._sff_segment.shapes:
                transform = self._transforms[shape.transformId]
                self._vtk_meshes.append(
                    VTKMesh.from_shape(shape, self.colour, self._vtk_args, transform)
                    )
        elif self.primary_descriptor == "threeDVolume":
            self._vtk_meshes.append(
                VTKMesh.from_volume(self._lattice, self.colour, self._vtk_args)
                )

    def __iter__(self):
        return iter(self._vtk_meshes)

    def __len__(self):
        return len(self._vtk_meshes)

    @property
    def primary_descriptor(self):
        return self._vtk_args.primary_descriptor

    @property
    def colour(self):
        return self._colour


class VTKContours(Contours):
    def __init__(self, mesh):
        self._mesh = mesh
        self._vtk_args = self._mesh.vtk_args
        self._colour = self._mesh.colour
        self._alpha = self._mesh.alpha
        self._mesh.ComputeBounds
        self.Xmin, self.Xmax, self.Ymin, self.Ymax, self.Zmin, self.Zmax = self._mesh.GetBounds()
        # create orthogonal planes
        self._axes = ['x', 'y', 'z']
        for d in self._axes:
            self.__setattr__('{}_plane'.format(d), vtk.vtkPlane())
            self.__getattribute__('{}_plane'.format(d)).SetOrigin(
                math.floor(self.Xmin),
                math.floor(self.Ymin),
                math.floor(self.Zmin),
                )
            d_hat = [0, 0, 0]
            d_hat[self._axes.index(d)] = 1
            self.__getattribute__('{}_plane'.format(d)).SetNormal(*d_hat)
        self._vtkcontours = self._build_vtkcontours()
        self._contours = self._build_contours()
        self._x_vtkcontours = None
        self._y_vtkcontours = None
        self._z_vtkcontours = None
        self._x_contours = None
        self._y_contours = None
        self._z_contours = None

    def _build_vtkcontours(self):
        # fill  holes
        filledHoles = vtk.vtkFillHolesFilter()
        filledHoles.SetInputData(self._mesh)
        filledHoles.SetHoleSize(1000)
        filledHoles.Update()
        filledHolesOutput = filledHoles.GetOutput()
        # first generate normals
        self._normals = vtk.vtkPolyDataNormals()
        self._normals.SetInputData(filledHolesOutput)
        self._normals.Update()
        cut_meshes = dict()
        for d in self._axes:
            """
            :TODO: parallelise
            """
            # the configure the cutter
            self._cutter = vtk.vtkCutter()
            self._cutter.SetInputConnection(self._normals.GetOutputPort())
            self._cutter.SetCutFunction(self.__getattribute__('%s_plane' % d))
            # determine the number of z-contours to generate from the bounds
            d_max = math.floor(self.__getattribute__('%smax' % d.upper()))
            d_count = int(2 * d_max + 1)
            # generate the cuts on the mesh
            self._cutter.GenerateValues(d_count, -d_max, d_max)
            # pass through vtk.vtkStripper to ensure we have triangel strips
            self._stripper = vtk.vtkStripper()
            self._stripper.SetInputConnection(self._cutter.GetOutputPort())
            self._stripper.Update()
            # the actual cuts
            colour = self._mesh.colour
            args = self._mesh.vtk_args
            cut_meshes[d] = VTKMesh(colour, args)
            stripper_output = self._stripper.GetOutput()
            cut_meshes[d].SetPoints(stripper_output.GetPoints())
            cut_meshes[d].SetPolys(stripper_output.GetLines())
        return cut_meshes

    def _build_contours(self):
        contour_dict = dict()
        for d in self._axes:
            contours = list()
            for i in xrange(self._vtkcontours[d].GetNumberOfCells()):
                cell = self._vtkcontours[d].GetCell(i)
                idList = cell.GetPointIds()
                contour = dict()
                for j in xrange(idList.GetNumberOfIds()):
                    contour[j] = self._vtkcontours[d].GetPoint(idList.GetId(j))
                contours.append(contour)
            contour_dict[d] = contours
        return contour_dict

    @property
    def colour(self):
        return self._colour

    @property
    def alpha(self):
        return self._alpha

    @property
    def x_vtkcontours(self):
        return self._vtkcontours['x']

    @property
    def y_vtkcontours(self):
        return self._vtkcontours['y']

    @property
    def z_vtkcontours(self):
        return self._vtkcontours['z']

    @property
    def x_contours(self):
        return self._contours['x']

    @property
    def y_contours(self):
        return self._contours['y']

    @property
    def z_contours(self):
        return self._contours['z']

    def render(self, renderer):
        assert isinstance(renderer, vtk.vtkRenderer)
        if any([
            self._vtk_args.all_contours,
            self._vtk_args.x_contours,
            self._vtk_args.y_contours,
            self._vtk_args.z_contours
            ]):
            _d = list()
            if self._vtk_args.all_contours:
                _d = ['x', 'y', 'z']
            else:
                if self._vtk_args.x_contours:
                    _d += ['x']
                if self._vtk_args.y_contours:
                    _d += ['y']
                if self._vtk_args.z_contours:
                    _d += ['z']
        for d in _d:
            self.mapper = vtk.vtkOpenGLPolyDataMapper()
            self.mapper.SetInputData(getattr(self, '{}_vtkcontours'.format(d)))
            self.actor = vtk.vtkOpenGLActor()
            self.actor.SetMapper(self.mapper)
            self.actor.GetProperty().SetColor(1, 1, 0)
            self.actor.GetProperty().SetOpacity(1)
            renderer.AddActor(self.actor)
        return renderer


class VTKSegment(Segment):
    def __init__(self, sff_segment, args, *args_, **kwargs_):
        self._sff_segment = sff_segment
        self._vtk_args = args
        self._meshes = VTKMeshes(self._sff_segment, self.colour, self._vtk_args, *args_, **kwargs_)
        self._contours = list()

    @property
    def id(self):
        return self._sff_segment.id

    @property
    def colour(self):
        """
        :TODO: currently only handles only RGBA colours
        """
        colour = list(self._sff_segment.colour.value)
        if colour[0] is not None and colour[1] is not None and colour[2] is not None and colour[3] is not None:
            colour[3] = self._vtk_args.transparency
            return colour
        else:
            colour = random(), random(), random(), 1
            print_date("Warning: random colour {} for segment {}".format(tuple(map(lambda x: round(x, 4), colour)), self._sff_segment.id))
            return colour

    @property
    def meshes(self):
        return self._meshes

    @property
    def contours(self):
        return self._contours

    def slice(self):
        for mesh in self.meshes:
            self._contours.append(mesh.slice())

    def render(self, renderer):
        for mesh in self.meshes:
            renderer = mesh.render(renderer)
        if any([self._vtk_args.all_contours, self._vtk_args.x_contours, self._vtk_args.y_contours, self._vtk_args.z_contours]):
            for contour_set in self.contours:
                renderer = contour_set.render(renderer)
        return renderer


class VTKHeader(Header):
    def __init__(self, sff_segment, args):
        self._sff_segment = sff_segment
        self._vtk_args = args


class VTKSegmentation(Segmentation):
    def __init__(self, sff_seg, args, configs):
        self._sff_seg = sff_seg  # the EMDB-SFF segmentation
        if not args.primary_descriptor:
            args.primary_descriptor = self._sff_seg.primaryDescriptor
        self._vtk_args = args
        self.configs = configs
        self._header = VTKHeader(self._sff_seg, self._vtk_args)
        self._segments = list()
        self._sliced_segments = list()
        if self._vtk_args.primary_descriptor == "threeDVolume":
            self._lattices = dict()
            for lattice in self._sff_seg.lattices:
                self._lattices[lattice.id] = lattice
                self._lattices[lattice.id].decode()
            for segment in self._sff_seg.segments:
                lattice = self._lattices[segment.volume.latticeId]
                lattice_data = lattice.data
                # new mask
                new_simplified_mask = numpy.ndarray(lattice_data.shape, dtype=int)
                new_simplified_mask[:,:,:] = 0
                # only the parts for this segment
                new_simplified_mask = (lattice_data == int(segment.volume.value)) * int(segment.volume.value)
                self._segments.append(
                    VTKSegment(
                        segment, self._vtk_args,
                        transforms=self._sff_seg.transforms,
                        lattice=new_simplified_mask
                    )
                )

            # if this is segmentation has threeDVolumes...
            # we need to build a mask for each segment
            # volume = self._sff_seg.segments[0].volume
            # with h5py.File(os.path.join(self._sff_seg.filePath, volume.file)) as f:
            #     try:
            #         # segger HDF5
            #         r_ids = f['/region_ids'].value  # placed at the top to fail fast if an EMDB map
            #         p_ids = f['/parent_ids'].value
            #         mask = f['/mask'].value
            #         c_ids = f['/region_colors'].value
            #         r_p_zip = zip(r_ids, p_ids)
            #         r_c_zip = dict(zip(r_ids, c_ids))
            #         simplified_mask, segment_ids = simplify_mask(mask, r_ids, r_p_zip, r_c_zip)
            #         for segment in self._sff_seg.segments:
            #             if segment.id in segment_ids:
            #                 new_simplified_mask = numpy.ndarray(mask.shape, dtype=int)
            #                 new_simplified_mask = 0
            #                 new_simplified_mask = (simplified_mask == segment.id) * segment.id
            #                 segment.mask = new_simplified_mask  # set the empty 'mask' attribute
            #                 self._segments.append(VTKSegment(segment, self._vtk_args))
            #     except KeyError:  # EMDB map HDF5
            #         for segment in self._sff_seg.segments:
            #             segment.mask = f['/mask'].value
            #             self._segments.append(VTKSegment(segment, self._vtk_args))
        else:
            self._segments = map(lambda s: VTKSegment(s, self._vtk_args, transforms=self._sff_seg.transforms), self._sff_seg.segments)

    @property
    def vtk_args(self):
        return self._vtk_args

    @property
    def header(self):
        return self._header

    @property
    def segments(self):
        return self._segments

    @property
    def bounds(self):
        _bounds = {
            'xmins': [],
            'xmaxs': [],
            'ymins': [],
            'ymaxs': [],
            'zmins': [],
            'zmaxs': []
            }
        for segment in self.segments:
            for mesh in segment.meshes:
                _xmin, _xmax, _ymin, _ymax, _zmin, _zmax = mesh.GetBounds()
                _bounds['xmins'].append(_xmin)
                _bounds['xmaxs'].append(_xmax)
                _bounds['ymins'].append(_ymin)
                _bounds['ymaxs'].append(_ymin)
                _bounds['zmins'].append(_zmin)
                _bounds['zmaxs'].append(_zmin)
        true_bounds = (
            min(_bounds['xmins']), max(_bounds['xmaxs']),
            min(_bounds['ymins']), max(_bounds['ymaxs']),
            min(_bounds['zmins']), max(_bounds['zmaxs'])
            )
        return true_bounds

    @property
    def center_point(self):
        xmin, xmax, ymin, ymax, zmin, zmax = self.bounds
        center = xmin + (xmax - xmin) / 2, ymin + (ymax - ymin) / 2, zmin + (zmax - zmin) / 2
        return center

    def slice(self):
        for segment in self.segments:
            segment.slice()

    def as_roi(self, configs):
        from ..formats.roi import ROISegmentation
        return ROISegmentation.from_vtk(self, configs)

    def render(self):
        """Render to display"""
        #  define the renderer
        ren = vtk.vtkOpenGLRenderer()
        ren.SetBackground(*self._vtk_args.background_colour)
        # populate the renderer with the meshes
        for segment in self.segments:
            ren = segment.render(ren)
        # render window
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        if self._vtk_args.full_screen:
            renWin.FullScreenOn()
        # render window interactor
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        # from vtkCubeAxesActor.h
        # define VTK_FLY_OUTER_EDGES     0
        # define VTK_FLY_CLOSEST_TRIAD   1
        # define VTK_FLY_FURTHEST_TRIAD  2
        # define VTK_FLY_STATIC_TRIAD    3
        # define VTK_FLY_STATIC_EDGES    4
        if self._vtk_args.cube_axes is not None:
            cubeAxesActor = vtk.vtkCubeAxesActor()
            cubeAxesActor.SetBounds(ren.ComputeVisiblePropBounds())
            cubeAxesActor.SetCamera(ren.GetActiveCamera())
            cubeAxesActor.SetFlyMode(self._vtk_args.cube_axes)
            cubeAxesActor.SetFlyModeToStaticEdges()  # how the cube axes will appear
            cubeAxesActor.GetTitleTextProperty(0).SetColor(1.0, 1.0, 1.0)
            cubeAxesActor.GetTitleTextProperty(1).SetColor(1.0, 1.0, 1.0)
            cubeAxesActor.GetTitleTextProperty(2).SetColor(1.0, 1.0, 1.0)
            cubeAxesActor.XAxisMinorTickVisibilityOff()
            cubeAxesActor.YAxisMinorTickVisibilityOff()
            cubeAxesActor.ZAxisMinorTickVisibilityOff()
            ren.AddActor(cubeAxesActor)
#             _actor_count += 1
#             assert ren.VisibleActorCount() == _actor_count
        # axes: display axes by default
        if not self._vtk_args.no_orientation_axes:
            axesActor = vtk.vtkAxesActor()
            axesWidget = vtk.vtkOrientationMarkerWidget()
            axesWidget.SetOrientationMarker(axesActor)
            axesWidget.SetViewport(0, 0, 0.1, 0.1)
            axesWidget.SetInteractor(iren)
            axesWidget.SetEnabled(1)
            ren.ResetCamera()
        # hello...
        print_date("Initialising...")
        iren.Initialize()
        print_date("Launching VTK viewer...")
        iren.Start()
        print_date("3D view completed.")

    def export(self, fn, args):
        json_data = dict()
        json_data['segments'] = list()
        center_point = self.center_point
        for segment in self.segments:
            segment_data = dict()
            segment_data['id'] = segment.id
            segment_data['colour'] = map(float, segment.colour)
            segment_data['meshes'] = list()
            for j, mesh in enumerate(segment.meshes):
                if args.center:
                    # center
                    x, y, z = center_point
                    to = (-x, -y, -z)
                    mesh2 = mesh.translate(to)
                else:
                    mesh2 = mesh
                # decimate
                decimateMesh = vtk.vtkDecimatePro()
                decimateMesh.SetInputData(mesh2)
                decimateMesh.SetTargetReduction(0.9)
                decimateMesh.PreserveTopologyOn()
                decimateMesh.Update()
                decimateMeshOutput = decimateMesh.GetOutput()
                writer = vtk.vtkXMLPolyDataWriter()
                out_fn = os.path.join(args.output_path, '{}_{}_m{}.vtp'.format(fn, segment.id, j))
                writer.SetFileName(out_fn)
                print_date("Exporting segment to {}".format(out_fn))
                writer.SetInputData(decimateMeshOutput)
                writer.SetDataModeToBinary()
                writer.SetHeaderTypeToUInt64()
                writer.SetCompressorTypeToZLib()
                writer.Write()
                segment_data['meshes'].append(os.path.basename(out_fn))
            json_data['segments'].append(segment_data)
        json_data['segment_count'] = len(self.segments)
        import json
        json_f = os.path.join(args.output_path, '{}_vtp_segments.json'.format(fn))
        with open(json_f, 'w') as f:
            json.dump(json_data, f, indent=4, sort_keys=True)
            if args.verbose:
                print_date("Exported metadata to {}".format(json_f))
