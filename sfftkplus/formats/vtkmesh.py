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

import json
import math
import os
from random import random

import numpy
import vtk
from sfftk.readers.segreader import get_root
from sfftkrw.core.print_tools import print_date, print_static
import sfftkrw as sff

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


class VTKMesh:
    def __init__(self, colour, args, *args_, **kwargs_):
        self._vtk_obj = vtk.vtkPolyData(*args_)
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
    def vtk_obj(self):
        return self._vtk_obj

    @vtk_obj.setter
    def vtk_obj(self, vtk_obj):
        try:
            assert isinstance(vtk_obj, vtk.vtkPolyData)
        except AssertionError:
            raise ValueError('vtk_obj must be of type vtk.vtkPolyData')
        self._vtk_obj = vtk_obj

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

    @staticmethod
    def transform_point(point, transform):
        """Apply the 3x4 transformation matrix to the 3-point

        :param point: tuple of 3 numbers
        :param transform: a 3x4 matrix
        :return point_transformed: the transformed point
        """
        point_transformed = transform.dot(numpy.array([list(point) + [1.0]]).T)
        return tuple(point_transformed.reshape(1, 3).flatten().tolist())

    @classmethod
    def from_mesh(cls, sff_mesh, colour, args, *args_, **kwargs_):
        """Initialiase a VTKMesh object from a ``sfftkrw.SFFMesh``

        :param mesh: a mesh with vertices and polygons
        :type mesh: ``sfftkrw.SFFMesh``
        :param colour: the segment colour
        :type colour: ``sfftkrw.SFFRGBA``
        :param args: parsed arguments
        :type args: ``argparse.Namespace``
        :return vtkmesh: an VTKMesh object
        :rtype vtkmesh: ``VTKMesh``
        """
        vertices = sff_mesh.vertices.data_array
        if sff_mesh.normals is not None:
            normals = sff_mesh.normals.data_array
        triangles = sff_mesh.triangles.data_array
        vtkmesh = cls(colour, args, *args_, **kwargs_)
        # define the geometry
        points = vtk.vtkPoints()
        for vertex_id, vertex in enumerate(vertices):
            points.InsertPoint(vertex_id, *vertex)
        vtkmesh.vtk_obj.SetPoints(points)
        # define the topology
        cellArray = vtk.vtkCellArray()
        for triangle in triangles:
            cell_size = len(triangle)
            cellArray.InsertNextCell(cell_size, triangle)
        vtkmesh.vtk_obj.SetPolys(cellArray)
        if not args.normals_off and sff_mesh.normals is not None:
            if len(normals) == len(vertices):
                _normals = vtk.vtkFloatArray()
                _normals.SetNumberOfComponents(3)
                for normal_id, normal in enumerate(normals):
                    _normals.InsertTuple3(normal_id, *normal)
                vtkmesh.vtk_obj.GetPointData().SetNormals(_normals)
        return vtkmesh

    @classmethod
    def from_volume(cls, lattice, colour, args, *args_, **kwargs_):
        """Initialiase a VTKMesh object from a binarised ``numpy.ndarray``
        
        :param lattice: a mask (3D matrix)
        :type lattice: ``numpy.ndarray``
        :param colour: the segment colour
        :type colour: ``sfftkrw.SFFRGBA``
        :param args: parsed arguments
        :type args: ``argparse.Namespace``
        :return vtkmesh: an VTKMesh object
        :rtype vtkmesh: ``VTKMesh``  
        """
        vtkmesh = cls(colour, args, *args_, **kwargs_)
        # try and figure out the lattice value if not specified
        if args.mask_value is not None:
            mask_value = args.mask_value
        else:
            mask_values = list(filter(lambda x: x != 0, set(lattice.flatten().tolist())))
            if len(mask_values) == 1:
                mask_value = mask_values[0]
            else:
                # use the mean
                mask_value = sum(mask_values) / len(mask_values)

        vol = vtk.vtkImageData()
        z_size, y_size, x_size = lattice.shape
        vol.SetExtent(0, x_size - 1, 0, y_size - 1, 0, z_size - 1)
        # vtkmesh.vol_extent = 0, x_size - 1, 0, y_size - 1, 0, z_size - 1
        vol.SetOrigin(0.0, 0.0, 0.0)
        sp_x = 1.0  # /(x_size - 1)
        sp_y = 1.0  # /(y_size - 1)
        sp_z = 1.0  # /(z_size - 1)
        vol.SetSpacing(sp_x, sp_y, sp_z)
        vol.AllocateScalars(vtk.VTK_FLOAT, 1)
        voxels = vtk.vtkFloatArray()
        voxels.SetNumberOfComponents(1)
        for m in lattice.flatten().tolist():
            voxels.InsertNextValue(float(m))
        vol.GetPointData().SetScalars(voxels)
        # convert the volume to a surface
        contours = vtk.vtkContourFilter()
        contours.SetInputData(vol)
        contours.SetValue(0, mask_value)
        contours.Update()
        contoursOutput = contours.GetOutput()

        vtkmesh.vtk_obj.SetPoints(contoursOutput.GetPoints())
        vtkmesh.vtk_obj.SetPolys(contoursOutput.GetPolys())
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
            #             print(vtkmesh_n.GetPointData().GetNormals())
            vtkmesh.vtk_obj.SetPoints(normals.GetOutput().GetPoints())
            vtkmesh.vtk_obj.SetPolys(normals.GetOutput().GetPolys())
            vtkmesh.vtk_obj.GetPointData().SetNormals(normalsOutput.GetPointData().GetNormals())
            return vtkmesh
        else:
            vtkmesh.vtk_obj.SetPoints(triangleStripsOutput.GetPoints())
            vtkmesh.vtk_obj.SetPolys(triangleStripsOutput.GetPolys())
            return vtkmesh

    @classmethod
    def from_shape(cls, shape, colour, args, transform, resolution=20, *args_, **kwargs_):
        """Initialiase a VTKMesh object from a sfftkrw.SFFShape
        
        :param shapes: an iterable of shapes
        :type shapes: ``sfftkrw.SFFShapePrimitiveList
        :param colour: the segment colour
        :type colour: ``sfftkrw.SFFRGBA``
        :param args: parsed arguments
        :type args: ``argparse.Namespace``
        :param transform: transform bearing this shape's translation from the origin
        :type transform: ``sfftkrw.SFFTransformationMatrix``
        :param int resolution: mesh resolution
        :return vtkmesh: an VTKMesh object
        :rtype vtkmesh: ``VTKMesh``  
        """
        assert resolution > 0
        vtkmesh = cls(colour, args, *args_, **kwargs_)
        from sfftkrw import SFFEllipsoid, SFFCuboid, SFFCylinder, SFFCone, SFFSubtomogramAverage
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
        elif isinstance(shape, SFFSubtomogramAverage):
            print("nothing to do here...")
        T = transform.data_array
        vtk_shape.SetCenter(float(T[0, 3]), float(T[1, 3]), float(T[2, 3]))
        vtk_shape.Update()
        _vtkmesh = vtk_shape.GetOutput()
        # triangle filter
        triangleMesh = vtk.vtkTriangleFilter()
        triangleMesh.SetInputData(_vtkmesh)
        triangleMesh.Update()
        triangleMeshOutput = triangleMesh.GetOutput()
        vtkmesh.vtk_obj.SetPoints(triangleMeshOutput.GetPoints())
        vtkmesh.vtk_obj.SetPolys(triangleMeshOutput.GetPolys())
        return vtkmesh

    @classmethod
    def from_vtk_obj(cls, vtk_obj, colour, args, *args_, **kwargs_):
        """Initialise a VTKMesh object from a VTKPolyData object"""
        vtkmesh = cls(colour, args, *args_, **kwargs_)
        vtkmesh.vtk_obj = vtk_obj
        return vtkmesh

    def slice(self):
        return VTKContours(self)

    def translate(self, to):
        # the transform
        translation = vtk.vtkTransform()
        translation.Translate(*to)
        # transform filter
        transformFilter = vtk.vtkTransformFilter()
        transformFilter.SetInputData(self.vtk_obj)
        transformFilter.SetTransform(translation)
        transformFilter.Update()
        return VTKMesh.from_vtk_obj(transformFilter.GetOutput(), self.colour, self._vtk_args)

    def scale(self, by):
        scale = vtk.vtkTransform()
        scale.Scale(*by)
        transformFilter = vtk.vtkTransformFilter()
        transformFilter.SetInputData(self.vtk_obj)
        transformFilter.SetTransform(scale)
        transformFilter.Update()
        return VTKMesh.from_vtk_obj(transformFilter.GetOutput(), self.colour, self._vtk_args)

    def rotate(self, by):
        print(f"by = {by}")
        rotation_matrix = vtk.vtkMatrix4x4()
        rotation_matrix.DeepCopy(by)
        rotation = vtk.vtkTransform()
        rotation.SetMatrix(rotation_matrix)
        transformFilter = vtk.vtkTransformFilter()
        transformFilter.SetInputData(self.vtk_obj)
        transformFilter.SetTransform(rotation)
        transformFilter.Update()
        return VTKMesh.from_vtk_obj(transformFilter.GetOutput(), self.colour, self._vtk_args)

    def render(self, renderer):
        """
        Render the mesh
         
        :param renderer: a renderer into which an actor derived from this mesh will be added
        :type renderer: `vtk.vtkRenderer`
        """
        assert isinstance(renderer, vtk.vtkRenderer)
        actor = self.make_actor()
        renderer.AddActor(actor)
        return renderer

    def make_actor(self):
        self.mapper = vtk.vtkOpenGLPolyDataMapper()
        # normals
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(self.vtk_obj)
        normals.ComputePointNormalsOn()
        normals.SplittingOff()
        normals.AutoOrientNormalsOn()
        normals.ConsistencyOn()
        normals.Update()
        self.mapper.SetInputData(normals.GetOutput())
        self.actor = vtk.vtkOpenGLActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(self.colour)
        if self.vtk_args.wireframe:
            self.actor.GetProperty().SetRepresentationToWireframe()
        if self.alpha:
            self.actor.GetProperty().SetOpacity(self.alpha)
        if self.vtk_args.view_edges:
            self.actor.GetProperty().EdgeVisibilityOn()
            self.actor.GetProperty().SetEdgeColor(*self.edge_colour)
        # lighting and shading
        # self.actor.GetProperty().SetInterpolationToPhong()
        # self.actor.GetProperty().SetDiffuse(0.7)
        # self.actor.GetProperty().SetSpecular(0.5)
        # self.actor.GetProperty().SetSpecularPower(40)
        # self.actor.GetProperty().BackfaceCullingOn()
        # self.actor.GetProperty().FrontfaceCullingOn()
        return self.actor


class VTKMeshes(object):
    def __init__(self, sff_segment, colour, args, *args_, **kwargs_):
        self._sff_segment = sff_segment
        self._vtk_args = args
        self._vtk_meshes = list()
        self._colour = colour
        # transforms
        if 'transform_list' in kwargs_:
            self._transforms = kwargs_['transform_list']
        if 'lattice' in kwargs_:
            self._lattice = kwargs_['lattice']
        if self.primary_descriptor == "mesh_list":
            for mesh in self._sff_segment.mesh_list:
                self._vtk_meshes.append(
                    VTKMesh.from_mesh(mesh, self.colour, self._vtk_args, *args_, **kwargs_)
                )
        elif self.primary_descriptor == "shape_primitive_list":
            # todo: handle each shape type independently e.g. ellipsoid, subtomogram_average, etc.
            if self._sff_segment.shape_primitive_list.num_subtomogram_averages:
                # todo: compute the vector by which to translate the lattice to have its centroid at the origin
                origin_vector = numpy.array([
                    [self._lattice.start.value[0] + self._lattice.size.value[0] / 2],
                    [self._lattice.start.value[1] + self._lattice.size.value[1] / 2],
                    [self._lattice.start.value[2] + self._lattice.size.value[2] / 2],
                ])
                # todo: compute the isosurface of the lattice; this is the reference mesh
                reference_mesh = VTKMesh.from_volume(self._lattice.data_array, self.colour, self._vtk_args)
                for sta in self._sff_segment.shape_primitive_list:
                    if isinstance(sta, sff.SFFSubtomogramAverage):
                        transform = self._transforms.get_by_id(sta.transform_id)
                        # todo: compute the rotation matrix; this will be done for each transform in a for loop
                        rotation_matrix = numpy.zeros((4, 4))
                        rotation_matrix[:3, :3] = transform.data_array[:3, :3]
                        rotation_matrix[3, 3] = 1  # homogeneous coordinates
                        # todo: the final translation vector
                        translation_vector = transform.data_array[:3, 3] / 4 # fixme: need a param for binning factor
                        origin_mesh = reference_mesh.translate((-origin_vector).flatten().tolist())
                        rotated_mesh = origin_mesh.rotate(rotation_matrix.flatten().tolist())
                        positioned_mesh = rotated_mesh.translate(translation_vector.flatten().tolist())
                        print(f"positioned_mesh: {positioned_mesh, type(positioned_mesh)}")
                        bounds = [0.0] * 6
                        positioned_mesh.vtk_obj.GetCellsBounds(bounds)
                        print(f"bounds: {bounds}")
                        self._vtk_meshes.append(positioned_mesh)
            # todo: add a parallel if handle for other shapes i.e. ellipsoids, cuboids, cones, cylinders.
            if self._sff_segment.shape_primitive_list.num_ellipsoids:
                print("Not implemented yet! Please contact the dev team to integrate this functionality.")
            if self._sff_segment.shape_primitive_list.num_cuboids:
                print("Not implemented yet! Please contact the dev team to integrate this functionality.")
            if self._sff_segment.shape_primitive_list.num_cones:
                print("Not implemented yet! Please contact the dev team to integrate this functionality.")
            if self._sff_segment.shape_primitive_list.num_cylinders:
                print("Not implemented yet! Please contact the dev team to integrate this functionality.")
        elif self.primary_descriptor == "three_d_volume":
            self._vtk_meshes.append(
                VTKMesh.from_volume(self._lattice, self.colour, self._vtk_args)
            )

    def __iter__(self):
        return iter(self._vtk_meshes)

    def __len__(self):
        return len(self._vtk_meshes)

    def __getitem__(self, index):
        return self._vtk_meshes[index]

    @property
    def primary_descriptor(self):
        return self._vtk_args.primary_descriptor

    @property
    def colour(self):
        return self._colour


class VTKContours(object):
    def __init__(self, mesh):
        self._mesh = mesh
        self._vtk_args = self._mesh.vtk_args
        self._colour = self._mesh.colour
        self._alpha = self._mesh.alpha
        self._mesh.vtk_obj.ComputeBounds()
        self.Xmin, self.Xmax, self.Ymin, self.Ymax, self.Zmin, self.Zmax = self._mesh.vtk_obj.GetBounds()
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
        self._contour_set = self._build_contours()
        self._x_vtkcontours = None
        self._y_vtkcontours = None
        self._z_vtkcontours = None
        self._x_contours = None
        self._y_contours = None
        self._z_contours = None

    def _build_vtkcontours(self):
        # fill  holes
        filledHoles = vtk.vtkFillHolesFilter()
        filledHoles.SetInputData(self._mesh.vtk_obj)
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
            cut_meshes[d].vtk_obj.SetPoints(stripper_output.GetPoints())
            cut_meshes[d].vtk_obj.SetPolys(stripper_output.GetLines())
        return cut_meshes

    def _build_contours(self):
        contour_set = dict()
        # for each of x, y, and z
        for d in self._axes:
            # build a set of contours for each orientation
            contours = list()
            for i in range(self._vtkcontours[d].vtk_obj.GetNumberOfCells()):
                cell = self._vtkcontours[d].vtk_obj.GetCell(i)
                idList = cell.GetPointIds()
                contour = dict()
                for j in range(idList.GetNumberOfIds()):
                    contour[j] = self._vtkcontours[d].vtk_obj.GetPoint(idList.GetId(j))
                contours.append(contour)
            contour_set[d] = contours
        return contour_set

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
        return self._contour_set['x']

    @property
    def y_contours(self):
        return self._contour_set['y']

    @property
    def z_contours(self):
        return self._contour_set['z']

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


class VTKSegment(object):
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
            print_date("Warning: random colour {} for segment {}".format(tuple(map(lambda x: round(x, 4), colour)),
                                                                         self._sff_segment.id))
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
        if any([self._vtk_args.all_contours, self._vtk_args.x_contours, self._vtk_args.y_contours,
                self._vtk_args.z_contours]):
            for contour_set in self.contours:
                renderer = contour_set.render(renderer)
        return renderer

    def append_meshes(self, appender):
        for mesh in self.meshes:
            appender.AddInputData(mesh)
        return appender


class VTKHeader(object):
    def __init__(self, sff_segment, args):
        self._sff_segment = sff_segment
        self._vtk_args = args


def decode_lattice(lattice):
    """Decodes a lattice by calling the decode method on a lattice object

    Void function

    :param lattice: :py:class:``sfftk.schema.SFFLattice`` object
    """
    lattice.decode()
    return lattice


class VTKSegmentation(object):
    def __init__(self, sff_seg, args, configs):
        self._sff_seg = sff_seg  # the EMDB-SFF segmentation
        if not args.primary_descriptor:
            args.primary_descriptor = self._sff_seg.primary_descriptor
        self._vtk_args = args
        self.configs = configs
        self._header = VTKHeader(self._sff_seg, self._vtk_args)
        self._segments = list()
        self._sliced_segments = list()
        # 3D volume segmentations
        if self._vtk_args.primary_descriptor == "three_d_volume":
            self._lattices = dict()
            # reconstitute into a dict
            if args.verbose:
                print_static("Decoding lattices...")
            for lattice in self._sff_seg.lattices:
                if args.verbose:
                    print_static("Decoding lattice {}...".format(lattice.id))
                # lattice.decode()
                print_date('', incl_date=False)
                self._lattices[lattice.id] = lattice.data_array
            # now we have the lattices decoded into 3D volumes,
            # we need to compute the surfaces for each
            for segment in self._sff_seg.segments:
                lattice_data = self._lattices[segment.three_d_volume.lattice_id]
                voxel_values = set(lattice_data.flatten().tolist()).difference({0})
                if len(voxel_values) == 1:  # it's a binary mask
                    if args.verbose:
                        print_date("Binary lattice")
                    new_simplified_mask = lattice_data
                else:
                    # lattice_data = lattice.data
                    if args.verbose:
                        print_static("Non-binary lattice: segment label #{}".format(int(segment.three_d_volume.value)))
                    # new mask
                    new_simplified_mask = numpy.ndarray(lattice_data.shape, dtype=numpy.dtype(int))
                    # new_simplified_mask = lattice
                    new_simplified_mask[:, :, :] = 0
                    # only the parts for this segment
                    new_simplified_mask = (lattice_data == int(segment.three_d_volume.value)) * int(
                        segment.three_d_volume.value)
                    if new_simplified_mask.sum() == 0:
                        print_date('', incl_date=False)
                        print_date('No data found for segment {}'.format(segment.id))
                        continue
                self._segments.append(
                    VTKSegment(
                        segment, self._vtk_args,
                        transforms=self._sff_seg.transform_list,
                        lattice=new_simplified_mask
                    )
                )
            print_date('', incl_date=False)
        elif self._vtk_args.primary_descriptor == "shape_primitive_list":
            for segment in self._sff_seg.segment_list:
                self._segments.append(
                    VTKSegment(
                        segment, self._vtk_args,
                        transform_list=self._sff_seg.transform_list,
                        lattice=self._sff_seg.lattice_list[0]
                    )
                )
            # todo: compute the isosurface of the lattice
            # todo: check if there are any subtomogram averages in the shape list
            # segment = self._sff_seg.segments[0]
            # image_to_physical_transform = self._sff_seg.transform_list.get_by_id(0)
            # image_to_physical_transform.data_array = numpy.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
            # # get the first transform
            # transform = self._sff_seg.transform_list.get_by_id(1)
            # transform_without_translation = numpy.zeros((4, 4))
            # transform_without_translation[:3, :3] = transform.data_array[:3, :3]
            # transform_without_translation[3, 3] = 1
            # print(f"image_to_physical_transform: {image_to_physical_transform.data_array}")
            # colour = segment.colour.value
            # if segment.shape_primitive_list.num_subtomogram_averages:
            #     lattice = self._sff_seg.lattices[0]
            #     print(f"lattice: {lattice}")
            #     # fixme: if we are in image space ignore voxel sizes; otherwise use them
            #     vx, vy, vz = numpy.diag(image_to_physical_transform.data_array)
            #     print(f"voxel sizes: {vx}, {vy}, {vz}")
            #     print(f"lattice dimensions: {lattice.size.value}")
            #     print(f"lattice origin: {lattice.start.value}")
            #     translate_to_origin = numpy.array([
            #         [vx * (lattice.start.value[0] + lattice.size.value[0] / 2)],
            #         [vy * (lattice.start.value[1] + lattice.size.value[1] / 2)],
            #         [vz * (lattice.start.value[2] + lattice.size.value[2] / 2)]
            #     ])
            #     print(f"translate_to_origin: {(-translate_to_origin).flatten().tolist()}")
            #     mesh = VTKMesh.from_volume(lattice.data_array, colour, self._vtk_args)
            #     mesh_at_origin = mesh.translate((-translate_to_origin).flatten().tolist())
            #     print(f"mesh_at_origin: {mesh_at_origin.vtk_obj}")
            #     # rotate about origin
            #     print(f"augmented_transform: {transform_without_translation.flatten().tolist()}")
            #     mesh_rotated = mesh_at_origin.rotate(transform_without_translation.flatten().tolist())
            #     print(f"mesh_rotated: {mesh_rotated.vtk_obj}")
            #     translate_vector = transform.data_array[:3, 3].flatten().tolist()
            #     print(f"translate_vector: {translate_vector}")
            #     mesh_at_position = mesh_rotated.translate(translate_vector)
            #     bounds = [0.0, ] * 6
            #     print(f"mesh_at_position: {mesh_at_position.vtk_obj.GetCellsBounds(bounds)}")
            #     print(f"bounds: {bounds}")
                # for each subtomogram average, create a transformed mesh representing the
                # location of the lattice's isosurface
                # for sta in segment.shape_primitive_list:
                #     transform = self._sff_seg.transform_list.get_by_id(sta.transform_id)
                #     contour_level = sta.value
                #     translate_to_origin = numpy.array([[0, 0, 0]]) # fixme: wrong value for now
                #     translate_to = transform.data_array[:3, 3]
                #     rotation = numpy.zeros((4, 4))
                #     rotation[:3, :3] = transform.data_array[:3, :3]
                #     print(f"rotation: {rotation}")
                #     print(f"transform: {transform.data_array}")
                #     print(f"to: {translate_to}")
                #     print(f"sta: {sta}")
                #     print(f"contour_level: {contour_level}")
                #     self._vtk_args.mask_value = contour_level
                #     # first, translate the sta to have it center aligned with the origin
                #     # the dimensions of the origin are specified by the dimensions and location of the lattice (complex)
                #     # next, perform rotations only (3x3 upper left matrix)
                #     # finally, translate to the position indicated by the right most vector
                #     print(sta)
                #     mesh = VTKMesh.from_volume(lattice.data_array, colour, self._vtk_args)
                #     print(f"mesh: {mesh.vtk_obj}")
                #     moved_mesh = mesh.translate(translate_to)
                #     print(f"moved_mesh: {moved_mesh}")
                # self._segments.append(
                #     VTKSegment(
                #         segment, self._vtk_args,
                #         transforms=self._sff_seg.transform_list,
                #         lattice=lattice.data_array
                #     )
                # )
        else:
            self._segments = list(map(lambda s: VTKSegment(s, self._vtk_args, transforms=self._sff_seg.transform_list),
                                      self._sff_seg.segment_list))

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
                _xmin, _xmax, _ymin, _ymax, _zmin, _zmax = mesh.vtk_obj.GetBounds()
                _bounds['xmins'].append(_xmin)
                _bounds['xmaxs'].append(_xmax)
                _bounds['ymins'].append(_ymin)
                _bounds['ymaxs'].append(_ymin)
                _bounds['zmins'].append(_zmin)
                _bounds['zmaxs'].append(_zmin)
        try:
            true_bounds = (
                min(_bounds['xmins']), max(_bounds['xmaxs']),
                min(_bounds['ymins']), max(_bounds['ymaxs']),
                min(_bounds['zmins']), max(_bounds['zmaxs'])
            )
        except ValueError:
            true_bounds = None
        return true_bounds

    @property
    def center_point(self):
        center = 0.0, 0.0, 0.0
        if self.bounds:
            xmin, xmax, ymin, ymax, zmin, zmax = self.bounds
            center = xmin + (xmax - xmin) / 2, ymin + (ymax - ymin) / 2, zmin + (zmax - zmin) / 2
        return center

    def slice(self):
        for segment in self.segments:
            segment.slice()

    def as_roi(self, args, configs):
        from .roi import ROISegmentation
        return ROISegmentation.from_vtk(self, args, configs)

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

    def export(self, fn, args, configs):
        json_data = dict()
        json_data['segments'] = list()
        center_point = self.center_point
        for segment in self.segments:
            segment_data = dict()
            segment_data['id'] = segment.id
            segment_data['colour'] = list(map(float, segment.colour))
            segment_data['meshes'] = list()
            for j, mesh in enumerate(segment.meshes):
                if args.apply_original_transform:
                    transform = self._sff_seg.transform_list[0].data_array
                    # scale
                    s_x = transform[0, 0]
                    s_y = transform[1, 1]
                    s_z = transform[2, 2]
                    # translate
                    t_x = transform[0, 3]
                    t_y = transform[1, 3]
                    t_z = transform[2, 3]
                    mesh3 = VTKMesh(mesh.colour, mesh.vtk_args)
                    # must scale first
                    mesh3.vtk_obj = mesh.scale((s_x, s_y, s_z)).vtk_obj
                    # then translate
                    mesh2 = mesh3.translate((t_x, t_y, t_z)).vtk_obj
                elif args.center:
                    # center
                    x, y, z = center_point
                    to = (-x, -y, -z)
                    mesh2 = mesh.translate(to).vtk_obj
                else:  # no transform
                    mesh2 = mesh.vtk_obj
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
        json_f = os.path.join(args.output_path, '{}_vtp_segments.json'.format(fn))
        with open(json_f, 'w') as f:
            json.dump(json_data, f, indent=4, sort_keys=True)
            if args.verbose:
                print_date("Exported metadata to {}".format(json_f))
