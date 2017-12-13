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
import numpy
import os
from random import random

import h5py
import vtk

from sfftk.core.print_tools import print_date
from sfftk.readers.segreader import get_root
from sfftkplus.formats.base import Segmentation, Header, Segment, Contours, Mesh


__author__  = "Paul K. Korir, PhD"
__email__   = "pkorir@ebi.ac.uk, paul.korir@gmail.com"
__date__    = "2017-04-12"


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
    simplified_mask = numpy.ndarray(mask.shape, dtype=int)  # @UndefinedVariable @UnusedVariable
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
        # that will do hundreds of array-wide comparisons at a time.
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


class ContoursToSurface:
    """Converts a stack of contours into a capped surface.
    
    Copyright 2014 Boyeong Woo
    """      
    def __init__(self, mode="Manual", distance_factor=3.0, use_default=True, resolution0=5, resolution1=5):
        self.mode = mode
        self.distance_factor = distance_factor
        self.use_default = use_default
        self.resolution0 = resolution0
        self.resolution1 = resolution1
    # end __init__

    def run(self, input):
        """ Run the filter using the selected mode. 
        
        input is a vtkPolyData object
        """
        # Get the outputs.
        if self.mode == "Manual":
            output = self.runManual(input)
        elif self.mode == "PointWalk":
            output = self.runPointWalk(input)
        else: # Resample
            output = self.runResample(input)
                
        # Triangulate remnant polygons
        triangle = vtk.vtkTriangleFilter()  # @UndefinedVariable
        triangle.SetInputData(output)
        triangle.Update()
        output = triangle.GetOutput()
        
        # Clear temporary lines used for processing (keep only polygons)
        emptyLines = vtk.vtkCellArray()  # @UndefinedVariable
        output.SetLines(emptyLines)
        
        # return the surface
        return output
        # end run

    def runPointWalk(self, input):
        """ Pass the input through the vtkRuledSurfaceFilter using the PointWalk mode. """

        surfaceFilter = vtk.vtkRuledSurfaceFilter()  # @UndefinedVariable

        if vtk.vtkVersion.GetVTKMajorVersion() <= 5:  # @UndefinedVariable
            surfaceFilter.SetInput(input)
        else:
            surfaceFilter.SetInputData(input)
        surfaceFilter.PassLinesOn() # necessary
        surfaceFilter.SetRuledModeToPointWalk()

        # Set the distance factor. It is used to decide when two lines are too far apart to connect.
        # The default distance factor set by the filter is 3.
        surfaceFilter.SetDistanceFactor(self.distance_factor)

        # Invoking OrientLoopsOn() causes the slicer to crash. See the report.
        #surfaceFilter.OrientLoopsOn()

        surfaceFilter.Update()
        return surfaceFilter.GetOutput()
    # end runPointWalk

    def runResample(self, input):
        """ Pass the input through the vtkRuledSurfaceFilter using the Resample mode. """
        surfaceFilter = vtk.vtkRuledSurfaceFilter()  # @UndefinedVariable

        if vtk.vtkVersion.GetVTKMajorVersion() <= 5:  # @UndefinedVariable
            surfaceFilter.SetInput(input)
        else:
            surfaceFilter.SetInputData(input)
        surfaceFilter.PassLinesOn() # necessary
        surfaceFilter.SetRuledModeToResample()

        # Set the resolutions.
        # resolution[0] defines the resolution in the direction of the parallel lines.
        # resolution[1] defines the resolution across the parallel lines.
        if self.use_default:
            # Use the maximum number of points on a line as resolution[0].
            # Default resolution[1] is 1.
            resolution0 = 0
            for i in range(input.GetNumberOfLines()):
                n = input.GetCell(i).GetNumberOfPoints()
                if n > resolution0:
                    resolution0 = n
            surfaceFilter.SetResolution(resolution0, 1)
        else:
            # Use the values specified by the user.
            surfaceFilter.SetResolution(self.resolution0, self.resolution1)
        surfaceFilter.Update()

        # Resample mode generates triangle strips.
        # Pass through the vtkTriangleFilter.
        triangleFilter = vtk.vtkTriangleFilter()  # @UndefinedVariable

        if vtk.vtkVersion.GetVTKMajorVersion() <= 5:  # @UndefinedVariable
            triangleFilter.SetInput(surfaceFilter.GetOutput())
        else:
            triangleFilter.SetInputData(surfaceFilter.GetOutput())

        triangleFilter.Update()
        return triangleFilter.GetOutput()
    # end runResample

    def runManual(self, originalInput):
        """ Try the manual PointWalk algorithm instead of the vtkRuledSurfaceFilter. """

        input = vtk.vtkPolyData() # @UndefinedVariable
        input.DeepCopy(originalInput)
        
        points = input.GetPoints()
        lines = input.GetLines()
        # add triangles to this
        polys = vtk.vtkCellArray() # @UndefinedVariable 
        
        # a line in a plane
        line1 = vtk.vtkLine() # @UndefinedVariable
        # a line in the next plane
        line2 = vtk.vtkLine() # @UndefinedVariable

        nlines = input.GetNumberOfLines() # total number of lines
                
        # remove keyholes from the lines
        newLines = self.fixKeyholes(input, nlines, 0.1, 2)
                
        self.setLinesClockwise(input, newLines)
        

        lines = vtk.vtkCellArray() # @UndefinedVariable
        lines.Initialize()
        input.DeleteCells()
        for i in newLines:
            lines.InsertNextCell(i)
        input.SetLines(lines)
        input.BuildCells()

        nlines = input.GetNumberOfLines() # total number of lines

        # Get two consecutive planes.        
        p1 = 0 # pointer to first line on plane 1
        nlns1 = self.getNumLinesOnPlane(input, nlines, p1) # number of lines on plane 1

        while p1 + nlns1 < nlines:

            p2 = p1 + nlns1 # pointer to first line on plane 2
            nlns2 = self.getNumLinesOnPlane(input, nlines, p2) # number of lines on plane 2

            # Initialize overlaps lists. - list of list
            # Each internal list represents a line from the plane and will store the pointers to the overlapping lines.

            # overlaps for lines from plane 1
            overlaps1 = []
            for i in range(nlns1):
                overlaps1.append([])

            # overlaps for lines from plane 2
            overlaps2 = []
            for j in range(nlns2):
                overlaps2.append([])

            # Fill the overlaps lists.
            for i in range(nlns1):
                line1.DeepCopy(input.GetCell(p1+i))
                for j in range(nlns2):
                    line2.DeepCopy(input.GetCell(p2+j))
                    if self.overlap(line1, line2):
                        # line i from plane 1 overlaps with line j from plane 2
                        overlaps1[i].append(p2+j)
                        overlaps2[j].append(p1+i)

            # Go over the overlaps lists.
            for i in range(p1, p1+nlns1):
                line1.DeepCopy(input.GetCell(i))
                pts1 = line1.GetPointIds()
                npts1 = line1.GetNumberOfPoints()

                intersects = False

                for j in overlaps1[i-p1]: # lines on plane 2 that overlap with line i

                    line2.DeepCopy(input.GetCell(j))
                    pts2 = line2.GetPointIds()
                    npts2 = line2.GetNumberOfPoints()

                    # Deal with possible branches.
                    #TODO: Improve the method for branching.
                    # (It only works for simple branching, so try some more complicated ones.)

                    # Get the portion of line 1 that is close to line 2.
                    divided1 = self.branch(input, pts1, npts1, j, overlaps1[i-p1])
                    dpts1 = divided1.GetPointIds()
                    dnpts1 = divided1.GetNumberOfPoints()

                    # Get the portion of line 2 that is close to line 1.
                    divided2 = self.branch(input, pts2, npts2, i, overlaps2[j-p2])
                    dpts2 = divided2.GetPointIds()
                    dnpts2 = divided2.GetNumberOfPoints()

                    # Use the divided lines for triangulation.
                    if dnpts1 > 1 and dnpts2 > 1:
                        self.manualPointWalk(input, dpts1, dnpts1, dpts2, dnpts2, polys, points)

                # end for j

            #end for i

            # Advance the pointers.
            p1 = p2
            nlns1 = nlns2

        # Triangulate all contours which are exposed.
        self.sealMesh(input, lines, polys)
        
        # Initialize the output data.
        output = vtk.vtkPolyData() # @UndefinedVariable
        output.SetPoints(points)
        output.SetLines(lines)
        # Return the output data.
        output.SetPolys(polys)
        return output
    # end runManual

    def fixKeyholes(self, input, nlines, epsilon, minSeperation):
        """ Parameters:
            - input: the input poly data
            - nlines: total number of lines in the input
            - epsilon: the threshold distance for keyholes
            - minSeperation: the minumum number of points of seperation between points
                                             for keyhole consideration
                                             
        Returns an array of vtkLine that represent the current structure with keyholes
        isolated
        """

        fixedLines = []
        for i in range(nlines):
            fixedLines += self.fixKeyhole(input, i, epsilon, minSeperation)
        return fixedLines

    def fixKeyhole(self, input, lineIndex, epsilon, minSeperation):
        """ Parameters:
            - input: the input poly data
            - line: the index of the current line
            - epsilon: the threshold distance for keyholes
            - minSeperation: the minumum number of points of seperation between points
                                             for keyhole consideration
                                             
        Returns an array of vtkLine that represents the current line with keyholes
        isolated
        """

        originalLine = vtk.vtkLine() # @UndefinedVariable
        originalLine.DeepCopy(input.GetCell(lineIndex))

        pts = originalLine.GetPoints()
        npts = originalLine.GetNumberOfPoints()

        # If the value of flags[i] is None, the point is not part of a keyhole
        # If the value of flags[i] is an integer, it represents a point that is
        # close enough that it could be considered part of a keyhole.
        flags = [None]*(npts-1)

        for i in range(npts-1):
            p1 = pts.GetPoint(i)

            for j in range(i+1, npts-1):
            
                # Make sure the points are not too close together on the line index-wise
                pointsOfSeperation = min(j-i,npts-1-j+i)
                if pointsOfSeperation > minSeperation:

                    # if the points are close together, mark both of them as part of a keyhole
                    p2 = pts.GetPoint(j)
                    distance = vtk.vtkMath().Distance2BetweenPoints(p1,p2) # @UndefinedVariable
                    if distance <= epsilon:
                        flags[i] = j
                        flags[j] = i

        newLines = []
        rawNewPoints = []
        finishedNewPoints = []

        currentLayer = 0

        inChannel = False

        # Loop through all of the points in the line
        for i in range(len(flags)):

            # Add a new line if neccessary
            if currentLayer == len(rawNewPoints):
                newLine = vtk.vtkLine() # @UndefinedVariable

                newLinePoints = newLine.GetPointIds()
                newLinePoints.Initialize()

                points = newLine.GetPoints()
                points.SetData(pts.GetData())

                newLines.append(newLine)
                rawNewPoints.append(newLinePoints)

            # If the current point is not part of a keyhole, add it to the current line
            if flags[i] == None:
                rawNewPoints[currentLayer].InsertNextId(originalLine.GetPointId(i))
                inChannel = False

            else:
                # If the current point is the start of a keyhole add the point to the line,
                # increment the layer, and start the channel.
                if flags[i] > i and not inChannel:
                    rawNewPoints[currentLayer].InsertNextId(originalLine.GetPointId(i))
                    currentLayer += 1
                    inChannel = True

                # If the current point is the end of a volume in the keyhole, add the point
                # to the line, remove the current line from the working list, deincrement
                # the layer, add the current line to the finished lines and start the,
                # channel.
                elif flags[i] < i and not inChannel:
                    rawNewPoints[currentLayer].InsertNextId(originalLine.GetPointId(i))
                    finishedNewPoints.append(rawNewPoints.pop())
                    currentLayer -= 1
                    inChannel = True

        # Add the remaining line to the finished list.
        for i in rawNewPoints:
            finishedNewPoints.append(i)

        # Seal the lines.
        for i in finishedNewPoints:
            if not i.GetNumberOfIds() == 0:
                i.InsertNextId(i.GetId(0))
        
        return newLines

    def setLinesClockwise(self, input, lines):
        ''' Parameters:
        - input: the input poly data
        - lines: the list of lines to be set clockwise
        
        Alters the lines in the list so that their points are oriented clockwise.
        '''
        
        numberOfLines = len(lines)
        
        for lineIndex in range(numberOfLines):
            if not self.clockwise(input, lines[lineIndex]):
                lines[lineIndex] = self.reverseLine(input, lines[lineIndex])
                
    # end setLinesClockwise

    def clockwise(self, input, line):
        ''' Parameters:
            - input: the input poly data
            - line: the vtkLine that is to be checked
            
            Returns 'True' if the specified line is oriented in a clockwise direction.
            'False' otherwise. Based on the shoelace algorithm.
        '''

        numberOfPoints = line.GetNumberOfPoints()

        # Calculate twice the area of the contour
        sum = 0
        for pointIndex in range(numberOfPoints-1):
            point1 = input.GetPoint(line.GetPointId(pointIndex))
            point2 = input.GetPoint(line.GetPointId(pointIndex+1))
            sum += (point2[0]-point1[0])*(point2[1]+point1[1])

        # If the area is positive, the contour is clockwise,
        # if it is negative, the contour is counter-clockwise.
        return sum > 0
    # end clockwise
    
    def reverseLine(self, input, originalLine):
        ''' Parameters:
            - input: the input poly data
            - originalLine: the vtkLine that is to be reversed
            
            Returns the vtkLine that is the reverse of the original.
        '''
        
        numberOfPoints = originalLine.GetNumberOfPoints()
        
        newLine = vtk.vtkLine() # @UndefinedVariable
        newPoints = newLine.GetPointIds()
        newPoints.Initialize()
        
        for pointInLineIndex in list(reversed(range(numberOfPoints))):
            newPoints.InsertNextId(originalLine.GetPointId(pointInLineIndex))

        return newLine
    # end reverseLine
        
    # helper function for runManual
    def getNumLinesOnPlane(self, input, nlines, p):
        """ Parameters:
            - input: the input poly data
            - nlines: total number of lines in the input
            - p: pointer to the first line on the plane of interest
            
        Returns the number of lines on the plane of interest.
        """
        plane = input.GetCell(p).GetBounds()[4] # z-value
        # Advance the pointer until z-value changes.
        i = p+1
        while i < nlines and input.GetCell(i).GetBounds()[4] == plane:
            i += 1
        return i-p
    # end getNumLinesOnPlane

    # helper function for runManual
    def overlap(self, line1, line2):
        """ Parameters:
            - line1 and line2: the two lines of interest
            
        Returns true if the bounds of the two lines overlap.
        """
        bounds1 = line1.GetBounds()
        bounds2 = line2.GetBounds()
        # true if there are overlaps in x-value ranges and y-value ranges of the two lines
        # bounds[0] is min-x; bounds[1] is max-x; bounds[2] is min-y; bounds[3] is max-y
        return bounds1[0] < bounds2[1] and bounds1[1] > bounds2[0] and bounds1[2] < bounds2[3] and bounds1[3] > bounds2[2]
    # end overlap

    # helper function
    def inside(self, line1, line2):
        """ Parameters:
            - line1 and line2: the two lines of interest
            
        Returns true if line1 is inside line2.
        (This function is not used yet. I'm not sure how this could be used and whether this would be useful or not.)
        """
        bounds1 = line1.GetBounds()
        bounds2 = line2.GetBounds()
        # true if x-value ranges and y-value ranges of line1 are inside those of line2
        return bounds1[0] > bounds2[0] and bounds1[1] < bounds2[1] and bounds1[2] > bounds2[2] and bounds1[3] < bounds2[3]
    # end inside

    # helper function for runManual
    def branch(self, input, pts, npts, i, overlaps):
        """ Parameters:
            - input: the input poly data
            - pts and npts: the line to be divided (trunk)
            - i: pointer to the line that is to be connected to the trunk (branch)
            - overlaps: list of all the lines that overlap with the trunk (possible branches)
            
        Get the portion of the trunk closest to the branch of interest.
        """
        divided = vtk.vtkLine() # @UndefinedVariable
        dpts = divided.GetPointIds()
        dpts.Initialize()

        # Discard some points on the trunk so that the branch connects to only a part of the trunk.
        prev = False # whether or not the previous point was added
        for j in range(npts):
            point = input.GetPoint(pts.GetId(j))
            # See if the point's closest branch is the input branch.
            if self.getClosestBranch(input, point, overlaps) == i:
                dpts.InsertNextId(pts.GetId(j))
                prev = True
            else:
                if prev:
                    # Add one extra point to close up the surface.
                    # (I'm not sure if this is the best way. You can change this if you want.)
                    dpts.InsertNextId(pts.GetId(j))
                prev = False
        dnpts = divided.GetNumberOfPoints()

        if dnpts > 1:
            # Determine if the trunk was originally a closed contour.
            closed = (pts.GetId(0) == pts.GetId(npts-1))
            if closed and (dpts.GetId(0) != dpts.GetId(dnpts-1)):
                # Make the new one a closed contour as well.
                # (I'm not sure if we have to make it closed always.)
                dpts.InsertNextId(dpts.GetId(0))

        return divided
    # end branch

    # helper function for branch
    def getClosestBranch(self, input, point, overlaps):
        """ Parameters:
            - input: the input poly data
            - point: a point on the trunk
            - overlaps: list of all the lines that overlap with the trunk
            
        Returns the branch closest from the point on the trunk.
        """
        # a branch from the overlaps
        line = vtk.vtkLine() # @UndefinedVariable 
        best = 0 # minimum distance from the point to the closest branch
        closest = 0 # pointer to the closest branch
        for i in overlaps:
            line.DeepCopy(input.GetCell(i))
            pts = line.GetPointIds()
            npts = line.GetNumberOfPoints()

            x = input.GetPoint(pts.GetId(0)) # a point from the branch
            # minimum distance from the point to the branch
            minD2 = vtk.vtkMath().Distance2BetweenPoints(point,x) # @UndefinedVariable 
            for j in range(1, npts):
                x = input.GetPoint(pts.GetId(j))
                distance2 = vtk.vtkMath().Distance2BetweenPoints(point,x) # @UndefinedVariable
                if distance2 < minD2:
                    minD2 = distance2

            # See if this branch is closer than the current closest.
            if best == 0 or minD2 < best:
                best = minD2
                closest = i

        return closest
    # end getClosestBranch

    # "polys" will be updated by this function (passed by reference).
    def manualPointWalk(self, input, pts1, npts1, pts2, npts2, polys, points):
        """This is modified from vtkRuledSurfaceFilter::PointWalk.
        
        Parameters:
            - input: the input poly data
            - pts1 and npts1: one of the lines to be connected (line 1)
            - pts2 and npts2: one of the lines to be connected (line 2)
            - polys: at the end, polys should hold the surface comprised of triangles
            
        :param input: input data
        :type input: `vtkPolyData`
        :param pts1: points of line1
        :type pts1: `vtkLine`
        :param int npts1: number of points in `pts1`
        :param pts2: points of line2
        :type pts2: `vtkLine`
        :param int npts2: number of points in `pts2`
        :param polys: resulting surface
        :type polys: `vtkPolyData`
        :param points: points in `polys`
        :type points: `vtkPoints`  
        """

        # Pre-calculate and store the closest points.

        # closest from line 1 to line 2
        closest1 = []
        for i in range(npts1):
            point = input.GetPoint(pts1.GetId(i))
            closest1.append(self.getClosestPoint(input, point, pts2, npts2))

        # closest from line 2 to line 1
        closest2 = []
        for i in range(npts2):
            point = input.GetPoint(pts2.GetId(i))
            closest2.append(self.getClosestPoint(input, point, pts1, npts1))

        # Orient loops.
        # Use the 0th point on line 1 and the closest point on line 2.
        # (You might want to try different starting points for better results.)
        startLoop1 = 0
        startLoop2 = closest1[0]

        x = input.GetPoint(pts1.GetId(startLoop1)) # first point on line 1
        a = input.GetPoint(pts2.GetId(startLoop2)) # first point on line 2
        xa = vtk.vtkMath().Distance2BetweenPoints(x,a) # @UndefinedVariable

        # Determine the maximum edge length.
        # (This is just roughly following the scheme from the vtkRuledSurfaceFilter.
        # You might want to try different numbers, or maybe not use this at all.)
        distance2 = xa * self.distance_factor * self.distance_factor

        # Determine if the loops are closed.
        # A loop is closed if the first point is repeated as the last point.
        # (I did not yet see any example where there were open contours, but I'll leave this for now.)
        loop1Closed = (pts1.GetId(0) == pts1.GetId(npts1-1))
        loop2Closed = (pts2.GetId(0) == pts2.GetId(npts2-1))

        # Determine the ending points.
        endLoop1 = self.getEndLoop(startLoop1, npts1, loop1Closed)
        endLoop2 = self.getEndLoop(startLoop2, npts2, loop2Closed)

        # for backtracking
        left = -1
        up = 1

        # Initialize the DP table.
        # Rows represent line 1. Columns represent line 2.

        # Fill the first row.
        firstRow = [xa]
        backtrackRow = [0]
        loc2 = self.getNextLoc(startLoop2, npts2, loop2Closed)
        for j in range(1, npts2):
            p = input.GetPoint(pts2.GetId(loc2)) # current point on line 2
            # Use the distance between first point on line 1 and current point on line 2.
            xp = vtk.vtkMath().Distance2BetweenPoints(x,p) # @UndefinedVariable
            firstRow.append(firstRow[j-1]+xp)
            backtrackRow.append(left)
            loc2 = self.getNextLoc(loc2, npts2, loop2Closed)

        # Fill the first column.
        score = [firstRow] # 2D list
        backtrack = [backtrackRow] # 2D list
        loc1 = self.getNextLoc(startLoop1, npts1, loop1Closed)
        for i in range(1, npts1):
            p = input.GetPoint(pts1.GetId(loc1)) # current point on line 1
            # Use the distance between first point on line 2 and current point on line 1
            pa = vtk.vtkMath().Distance2BetweenPoints(p,a) # @UndefinedVariable
            score.append([score[i-1][0]+pa]+[0]*(npts2-1)) # appending another row
            backtrack.append([up]+[0]*(npts2-1)) # appending another row
            loc1 = self.getNextLoc(loc1, npts1, loop1Closed)

        # Fill the rest of the table.
        prev1 = startLoop1
        prev2 = startLoop2
        loc1 = self.getNextLoc(startLoop1, npts1, loop1Closed)
        loc2 = self.getNextLoc(startLoop2, npts2, loop2Closed)
        for i in range(1, npts1):
            x = input.GetPoint(pts1.GetId(loc1)) # current point on line 1
            for j in range(1, npts2):
                a = input.GetPoint(pts2.GetId(loc2)) # current point on line 2
                xa = vtk.vtkMath().Distance2BetweenPoints(x,a) # @UndefinedVariable

                # Use the pre-calculated closest point.
                # (This was not in the report. Sometimes it seemed like it had to take a longer span, so that's why I added this,
                # but I have not yet fully tested this. If this does not seem much better, you can remove this part.)
                if loc1 == closest2[prev2]: # if loc1 is the closest from prev2
                    score[i][j] = score[i][j-1]+xa
                    backtrack[i][j] = left
                elif loc2 == closest1[prev1]: # if loc2 is the closest from prev1
                    score[i][j] = score[i-1][j]+xa
                    backtrack[i][j] = up

                # score[i][j] = min(score[i][j-1], score[i-1][j]) + xa
                elif score[i][j-1] <= score[i-1][j]:
                    score[i][j] = score[i][j-1]+xa
                    backtrack[i][j] = left
                else:
                    score[i][j] = score[i-1][j]+xa
                    backtrack[i][j] = up

                # Advance the pointers.
                prev2 = loc2
                loc2 = self.getNextLoc(loc2, npts2, loop2Closed)
            prev1 = loc1
            loc1 = self.getNextLoc(loc1, npts1, loop1Closed)

        # Backtrack.
        loc1 = endLoop1
        loc2 = endLoop2
        while i > 0 or j > 0:
            x = input.GetPoint(pts1.GetId(loc1)) # current point on line 1
            a = input.GetPoint(pts2.GetId(loc2)) # current point on line 2
            xa = vtk.vtkMath().Distance2BetweenPoints(x,a) # @UndefinedVariable

            if backtrack[i][j] == left:
                prev2 = self.getPrevLoc(loc2, npts2, loop2Closed)
                b = input.GetPoint(pts2.GetId(prev2)) # previous point on line 2
                xb = vtk.vtkMath().Distance2BetweenPoints(x,b) # @UndefinedVariable
                # Insert triangle if the spans are not larger than the maximum distance.

                #if xa <= distance2 and xb <= distance2:\
                currentTriangle = [pts1.GetId(loc1), pts2.GetId(loc2), pts2.GetId(prev2), pts1.GetId(loc1)]
                polys.InsertNextCell(3)
                polys.InsertCellPoint(currentTriangle[0])
                polys.InsertCellPoint(currentTriangle[1])
                polys.InsertCellPoint(currentTriangle[2])

                # Advance the pointers (backwards).
                j -= 1
                loc2 = prev2

            else: # up
                prev1 = self.getPrevLoc(loc1, npts1, loop1Closed)
                y = input.GetPoint(pts1.GetId(prev1)) # previous point on line 1
                ya = vtk.vtkMath().Distance2BetweenPoints(y,a) # @UndefinedVariable
                # Insert triangle if the triangle does not go out of the bounds of the contours.
                # If checkSurface returns None, the triangle is fine, otherwise it returns a list of new triangles
                # to be added instead
                currentTriangle = [pts1.GetId(loc1), pts2.GetId(loc2), pts1.GetId(prev1), pts1.GetId(loc1)]
                polys.InsertNextCell(3)
                polys.InsertCellPoint(currentTriangle[0])
                polys.InsertCellPoint(currentTriangle[1])
                polys.InsertCellPoint(currentTriangle[2])

                # Advance the pointers (backwards).
                i -= 1
                loc1 = prev1
    #TODO: Deal with rapid changes and internal contours.
    # (I was thinking it might be useful to be able to detect when a triangle goes outside the contour lines so that it can be corrected.)
    # end manualPointWalk

    def sealMesh(self, input, lines, polys):
        '''
        Parameters:
            - input: the input poly data
            - lines: vtkCellArray representing all the contours in the mesh
            - polys: vtkCellArray representing all of the polys in the mesh
            
        This function seals all contours that do not have polygons connecting on both sides.
        '''

        numLines = lines.GetNumberOfCells()
        numPolys = polys.GetNumberOfCells()

        # Keep track of whether polygons connect to contours from above or below.
        polygonsToBelow = [False]*numLines
        polygonsToAbove = [False]*numLines

        # Loop through the lines.
        line = vtk.vtkLine() # @UndefinedVariable
        for i in range(numLines-1):            
            line.DeepCopy(input.GetCell(i))
            pts = line.GetPointIds()
            z = line.GetBounds()[4]

            # Loop through the polygons
            polygonIds = vtk.vtkIdList() # @UndefinedVariable
            polys.SetTraversalLocation(0)
            for j in range(numPolys):

                # If polygons connect to the current contour from above and below,
                # it doesnt need to be sealed.
                if polygonsToAbove[i] and polygonsToBelow[i]:
                    break

                # Get the Id and z coordinates of the polygon corners.
                polys.GetNextCell(polygonIds)
                p1Id = polygonIds.GetId(0)
                p2Id = polygonIds.GetId(1)
                p3Id = polygonIds.GetId(2)
                z1 = input.GetPoint(p1Id)[2]
                z2 = input.GetPoint(p2Id)[2]
                z3 = input.GetPoint(p3Id)[2]

                # Only check polygons which lie on the same Z.
                if z1 == z or z2 == z or z3 == z:

                    # Check to see if the corners of the polygon lie on the line.
                    p1OnLine = self.onLine(pts, p1Id)
                    p2OnLine = self.onLine(pts, p2Id)
                    p3OnLine = self.onLine(pts, p3Id)

                    # If any of the corners of the current polygon lies on the current line.
                    if p1OnLine or p2OnLine or p3OnLine:
                        if not p1OnLine:
                            zNotOnLine = z1
                        elif not p2OnLine:
                            zNotOnLine = z2
                        elif not p3OnLine:
                            zNotOnLine = z3

                        # Check to see if the current polygon connects to the contour above or below it.
                        if zNotOnLine > z:
                            polygonsToAbove[i] = True
                        else:
                            polygonsToBelow[i] = True

        # Seal all contours which are only connected on one side.
        for i in range(numLines):
            line.DeepCopy(input.GetCell(i))
            if not (polygonsToAbove[i] and polygonsToBelow[i]):
                self.triangulateLine(input, line, polys)
    # end sealMesh

    def onLine(self, lineIds, id):
        '''
        Parameters:
        - lineIds: vtkIdList of the points in a line
        - id: id that is being checked
        Returns true if 'id' is in 'lineIds'
        '''

        numberOfPoints = lineIds.GetNumberOfIds()
        for currentId in range(numberOfPoints):
            if lineIds.GetId(currentId) == id:
                return True
        return False    
    # end onLine

    def triangulateLine(self, input, line, polys):    
        '''
        Parameters:
        - input: the input poly data
        - line: vtkLine representing the line to be triangulated
        - polys: vtkCellArray containing current polygons
        This function adds new polygons to polys that represent a 2D
        triangulation of 'line' on the x-y plane.
        '''

        npts = line.GetNumberOfPoints()

        linePolyData = vtk.vtkPolyData() # @UndefinedVariable
        linePolyData.SetPoints(line.GetPoints())

        boundary = vtk.vtkPolyData() # @UndefinedVariable
        boundary.SetPoints(linePolyData.GetPoints())

        # Use vtkDelaunay2D to triangulate the line
        # and produce new polyons
        delaunay = vtk.vtkDelaunay2D() # @UndefinedVariable

        if vtk.vtkVersion.GetVTKMajorVersion() <= 5: # @UndefinedVariable
            delaunay.SetInput(linePolyData)
            delaunay.SetSource(boundary)
        else:
            delaunay.SetInputData(linePolyData)
            delaunay.SetSourceData(boundary)
        delaunay.Update()
        out = delaunay.GetOutput()
        newPolygons = out.GetPolys()        

        # Check each new polygon to see if it is inside the line.
        ids = vtk.vtkIdList()      # @UndefinedVariable
        for i in range(out.GetNumberOfPolys()):
            newPolygons.GetNextCell(ids)

            # Get the center of the polygon
            x = 0
            y = 0
            for i in range(3):
                point = input.GetPoint(line.GetPointId(ids.GetId(i)))
                x += point[0]
                y += point[1]
            x /= 3
            y /= 3

            # Check to see if the center of the polygon lies inside the line.
            if self.pointInsideLine(input, line, (x,y)):
                polys.InsertNextCell(3)
                polys.InsertCellPoint(line.GetPointId(ids.GetId(0)))
                polys.InsertCellPoint(line.GetPointId(ids.GetId(1)))
                polys.InsertCellPoint(line.GetPointId(ids.GetId(2)))
    # end triangulateLine

    def pointInsideLine(self, input, line, point):
        '''
        Parameters:
        - input: the input poly data
        - line: vtkLine representing the line that is being checked
        - point: an (x, y) tuple that represents the point that is being checked
        This function uses ray casting to determine if a point lies within the line
        '''

        # Create a ray that starts outside the polygon and goes to the point being checked
        bounds = line.GetBounds()    
        rayPoint1 = (bounds[0]-10, bounds[2]-10)
        rayPoint2 = point
        ray = (rayPoint1,rayPoint2)

        # Check all of the edges to see if they intersect
        numberOfIntersections = 0
        for i in range(line.GetNumberOfPoints()-1):
            edgePoint1 = input.GetPoint(line.GetPointId(i))
            edgePoint2 = input.GetPoint(line.GetPointId(i+1))
            edge = (edgePoint1,edgePoint2)

            if point == edgePoint1 or point == edgePoint2:
                return True

            if self.lineIntersection(ray, edge, False):
                numberOfIntersections += 1

        # If the number of intersections is odd, the point is inside
        return numberOfIntersections%2 == 1
    # end pointInsideLine

    def lineIntersection(self, line1, line2, countInfinite):
        '''
        Parameters:
        - input: the input poly data
        - line1: tuple of ((lx1,ly1),(lx2,ly2)) that identifies the first line segment
        - line2: tuple of ((lx1,ly1),(lx2,ly2)) that identifies the second line segment
        - countInfinite: boolean which dictates if lines that are coincident should be considered to intersect
        This function returns True if the two line segments intersect, False otherwise.
        '''

        # Get the bounding region for the intersection
        xmin = max(min(line1[0][0],line1[1][0]),min(line2[0][0],line2[1][0]))
        xmax = min(max(line1[0][0],line1[1][0]),max(line2[0][0],line2[1][0]))
        ymin = max(min(line1[0][1],line1[1][1]),min(line2[0][1],line2[1][1]))
        ymax = min(max(line1[0][1],line1[1][1]),max(line2[0][1],line2[1][1]))

        # If the two lines don't overlap, no intersection
        if xmin > xmax or ymin > ymax:
            return False

        slope1 = self.getSlope(line1)
        slope2 = self.getSlope(line2)

        intercept1 = self.getIntercept(line1, slope1)
        intercept2 = self.getIntercept(line2, slope2)

        # Lines either parallel or coincident
        if (slope1 == slope2):
            if slope1 == None:
                return countInfinite and line1[0][1] == line2[0][1]
            else:
                return countInfinite and intercept1 == intercept2

        # If one of the lines has no distance, check to see if it lies on the other line as a point
        if line1[0] == line1[1] and line2[0] == line2[1]:
            return line1[0] == line2[0]
        elif line1[0] == line1[1]:
            return self.pointOnLine(slope2, intercept2, line1[0])
        elif line2[0] == line2[1]:
            return self.pointOnLine(slope1, intercept1, line2[0])

        # Find the intersection location for the lines
        if slope1 == None:
            x = line1[0][0]
            y = slope2*x+intercept2
        elif slope2 == None:
            x = line2[0][0]
            y = slope1*x+intercept1
        else:
            x = (intercept2-intercept1)/(slope1-slope2)
            if slope1 == 0:
                y = intercept1
            elif slope2 == 0:
                y = intercept2
            else:
                y = slope1*x+intercept1

        # If the intersection lies within the segments, return True, otherwise return False.
        if x >= xmin and\
             x <= xmax and\
             y >= ymin and\
             y <= ymax:         
            return True
        else:
            return False
    # end lineIntersection

    def getSlope(self, line):
        '''
        Parameters:
        - input: the input poly data
        - line: tuple of ((lx1,ly1),(lx2,ly2)) that identifies the line
        This function returns the slope of the line, or None if the slope is infinite.
        '''
        xChange = line[0][0] - line[1][0]
        if xChange == 0:
            return None
        else:
            yChange = line[0][1] - line[1][1]
            return yChange/xChange
    # end getSlope

    def getIntercept(self, edge, slope):
        '''
        Parameters:
        - input: the input poly data
        - line: tuple of ((lx1,ly1),(lx2,ly2)) that identifies the line
        This function returns the y-intercept of the line, or None if the slope is infinite.
        '''
        if slope == None:
            return None
        else:
            return edge[0][1] - slope*edge[0][0] 
    # end getIntercept

    def pointOnLine(self, m, b, point):
        '''
        Parameters:
        - m: the slope of the line
        - b: the y-intercept of the line
        - point: an (x, y) tuple that represents the point that is being checked
        This function returns True if the point lies on the line, False otherwise.
        '''
        epsilon = 0.00001
        y = m*point[0] + b
        return abs(y - point[1]) < epsilon
    # end pointOnLine             

    # helper function for manualPointWalk
    def getClosestPoint(self, input, point, pts, npts):
        """ Parameters:
            - input: the input poly data
            - point: a point from a line
            - pts and npts: the other line that is to be connected
            
        Returns the point on the given line that is closest to the given point. """
        x = input.GetPoint(pts.GetId(0)) # point from the given line
        # minimum distance from the point to the line
        minD2 = vtk.vtkMath().Distance2BetweenPoints(point,x) # @UndefinedVariable 
        closest = 0 # pointer to the closest point
        for i in range(1, npts):
            x = input.GetPoint(pts.GetId(i))
            distance2 = vtk.vtkMath().Distance2BetweenPoints(point,x) # @UndefinedVariable
            if distance2 < minD2:
                minD2 = distance2
                closest = i
        return closest
    # end getClosestPoint

    # helper function for manualPointWalk
    def getEndLoop(self, startLoop, npts, loopClosed):
        """ Parameters:
            - startLoop: pointer to the starting point
            - npts: number of points on the loop
            - loopClosed: whether or not the loop is closed
            
        Returns the ending point for the loop. """
        if startLoop != 0:
            if loopClosed:
                return startLoop
            return startLoop-1
        # If startLoop was 0, then it doesn't matter whether or not the loop was closed.
        return npts-1
    # end getEndLoop

    # helper function for manualPointWalk
    def getNextLoc(self, loc, npts, loopClosed):
        """ Parameters:
            - loc: pointer to the current point
            - npts: number of points on the loop
            - loopClosed: whether or not the loop is closed
            
        Returns the next point on the loop. """
        if loc+1 == npts: # if the current location is the last point
            if loopClosed:
                # Skip the repeated point.
                return 1
            return 0
        return loc+1
    # end getNextLoc

    # helper function for manualPointWalk
    def getPrevLoc(self, loc, npts, loopClosed):
        """ Parameters:
            - loc: pointer to the current point
            - npts: number of points on the loop
            - loopClosed: whether or not the loop is closed
            
        Returns the previous point on the loop. """
        if loc-1 == -1: # if the current location is the first point
            if loopClosed:
                # Skip the repeated point.
                return npts-2
            return npts-1
        return loc-1


class VTKMesh(vtk.vtkPolyData):  # @UndefinedVariable
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
                if P.vertex_ids[0] in normal_ids: # if the first vertex in the polygon is a normal
                    polygons[P.PID] = tuple(P.vertex_ids[1::2])
                elif P.vertex_ids[0] in vertex_ids: # if the first vertex in the polygon is a vertex
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
        points = vtk.vtkPoints()  # @UndefinedVariable
        for vertex_id, vertex in new_vertices.iteritems():
            points.InsertPoint(vertex_id, *vertex)
        vtkmesh.SetPoints(points)
        # define the topology
        cellArray = vtk.vtkCellArray()  # @UndefinedVariable
        for polygon in new_polygons.itervalues():
            cell_size = len(polygon)
            cellArray.InsertNextCell(cell_size, polygon)
        vtkmesh.SetPolys(cellArray)
        # define normals (if they exist)
        if not args.normals_off:
            if len(new_normals) == len(new_vertices):
                normals = vtk.vtkFloatArray()  # @UndefinedVariable
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
                mask_value = sum(mask_values)/len(mask_values)
                
        vol = vtk.vtkImageData()  # @UndefinedVariable
        z_size, y_size, x_size = mask.shape
        vol.SetExtent(0, x_size - 1, 0, y_size - 1 , 0, z_size - 1)
        vtkmesh.vol_extent = 0, x_size - 1, 0, y_size - 1 , 0, z_size - 1
        vol.SetOrigin(0.0, 0.0, 0.0)
        sp_x = 1.0 #/(x_size - 1)
        sp_y = 1.0 #/(y_size - 1)
        sp_z = 1.0 #/(z_size - 1)
        vol.SetSpacing(sp_x, sp_y, sp_z)
        vol.AllocateScalars(vtk.VTK_FLOAT, 1)  # @UndefinedVariable
        voxels = vtk.vtkFloatArray()  # @UndefinedVariable
        voxels.SetNumberOfComponents(1)
        for m in mask.flatten().tolist():
            voxels.InsertNextValue(float(m))
        vol.GetPointData().SetScalars(voxels)
        # convert the volume to a surface
        contours = vtk.vtkContourFilter()  # @UndefinedVariable
        contours.SetInputData(vol)
        contours.SetValue(0, mask_value)
        contours.Update()
        contoursOutput = contours.GetOutput()
        
        vtkmesh.SetPoints(contoursOutput.GetPoints())
        vtkmesh.SetPolys(contoursOutput.GetPolys())
        # triangulate
        triangleMesh = vtk.vtkTriangleFilter()  # @UndefinedVariable
        triangleMesh.SetInputData(contoursOutput)
        triangleMesh.Update()
        triangleMeshOutput = triangleMesh.GetOutput()
        # decimate
        decimateMesh = vtk.vtkDecimatePro()  # @UndefinedVariable
        decimateMesh.SetInputData(triangleMeshOutput)
        decimateMesh.SetTargetReduction(0.9)
        decimateMesh.PreserveTopologyOn()
        decimateMesh.Update()
        decimateMeshOutput = decimateMesh.GetOutput()
        # smooth
        smoothedMesh = vtk.vtkSmoothPolyDataFilter()  # @UndefinedVariable
        smoothedMesh.SetInputData(decimateMeshOutput)
        smoothedMesh.SetNumberOfIterations(200)
        smoothedMesh.Update()
        smoothedMeshOutput = smoothedMesh.GetOutput()
        # finally set things
        if not args.normals_off:
            # normals
            normals = vtk.vtkPolyDataNormals()
            normals.SetInputData(smoothedMeshOutput)
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
            vtkmesh.SetPoints(smoothedMeshOutput.GetPoints())
            vtkmesh.SetPolys(smoothedMeshOutput.GetPolys())
            return vtkmesh
    @classmethod
    def from_contours(cls, contours, colour, args):
        """Initialiase a VTKMesh object from a ``sfftk.schema.SFFContourList``
        
        :param contours: an iterable of contours
        :type contours: ``sfftk.schema.SFFContourList``
        :param colour: the segment colour
        :type colour: ``sfftk.schema.SFFColour``
        :param args: parsed arguments
        :type args: ``argparse.Namespace``
        :return vtkmesh: an VTKMesh object
        :rtype vtkmesh: ``VTKMesh``  
        """
        vtkmesh = cls(colour, args)
        # pack the contours as vertices and polygons
        vertex_id = 0
        vertices = dict()
        polygons = list()
        for contour in contours:
            polygon = list()
            for point in contour.points:
                vertices[vertex_id] = point.value
                polygon.append(vertex_id)
                vertex_id += 1
            polygons.append(polygon)
        # PROTOCOL
        # create a vtkAppendPolyData object to hold all the vtkPolyData objects
        # build each contour as a vtkPolyData and append
        # assign points and polys to self
        appendFilter = vtk.vtkAppendPolyData()  # @UndefinedVariable
        for polygon in polygons:
            if not polygon: # sometimes we have empty cell
                continue
            points = vtk.vtkPoints()  # @UndefinedVariable
            cells = vtk.vtkCellArray()  # @UndefinedVariable
            # init
            points.SetNumberOfPoints(len(polygon) + 1) # +1 because we close the cell
            cells.Allocate(1, len(polygon) + 1)
            cells.InsertNextCell(len(polygon) + 1)
            # add points
            for i in xrange(len(polygon)):
                points.SetPoint(i, *vertices[polygon[i]])
                cells.InsertCellPoint(i)
            # close the cells
            points.SetPoint(i + 1, *vertices[polygon[i]])
            cells.InsertCellPoint(i + 1)
            # make a polydata object from these primitives
            poly = vtk.vtkPolyData()  # @UndefinedVariable
            poly.Initialize()
            poly.SetPolys(cells)
            poly.SetPoints(points)
            # append
            appendFilter.AddInputData(poly)
        appendFilter.Update()
        appendOutput = appendFilter.GetOutput()
        # clean
        cleaned = vtk.vtkCleanPolyData()  # @UndefinedVariable
        cleaned.SetInputData(appendOutput)
        cleaned.Update() 
        cleanedOutput = cleaned.GetOutput()
        # depthsort
        depthsort = vtk.vtkDepthSortPolyData()  # @UndefinedVariable
        depthsort.SetInputData(cleanedOutput)
        depthsort.SetVector(0, 0, 1)
        depthsort.SetOrigin(0, 0, 0)
        depthsort.SetDirectionToSpecifiedVector()
        depthsort.SetDepthSortModeToParametricCenter()
        depthsort.Update()
        depthsortOutput = depthsort.GetOutput()
        # set primitives
        vtkmesh.SetPoints(depthsortOutput.GetPoints())
        vtkmesh.SetLines(depthsortOutput.GetPolys())
        # convert object
        C2S = ContoursToSurface()
        mesh = C2S.run(vtkmesh)
        # set primitives
        vtkmesh.SetPoints(mesh.GetPoints())
        vtkmesh.SetPolys(mesh.GetPolys())
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
            vtk_shape = vtk.vtkSphereSource()  # @UndefinedVariable
            vtk_shape.SetRadius(shape.x)
            """
            :TODO: make this variable
            """
            vtk_shape.SetPhiResolution(resolution)
            vtk_shape.SetThetaResolution(resolution)
        elif isinstance(shape, SFFCylinder):
            vtk_shape = vtk.vtkCylinderSource()  # @UndefinedVariable
            vtk_shape.SetHeight(shape.height)
            vtk_shape.SetRadius(shape.diameter/2)
            vtk_shape.SetResolution(resolution)
        elif isinstance(shape, SFFCone):
            vtk_shape = vtk.vtkConeSource()  # @UndefinedVariable
            vtk_shape.SetHeight(shape.height)
            vtk_shape.SetRadius(shape.bottomRadius)
            vtk_shape.SetResolution(resolution)
        elif isinstance(shape, SFFCuboid):
            vtk_shape = vtk.vtkCubeSource()  # @UndefinedVariable
            vtk_shape.SetXLength(shape.x)
            vtk_shape.SetYLength(shape.y)
            vtk_shape.SetZLength(shape.z)
        T = transform.data_array
        vtk_shape.SetCenter(float(T[0,3]), float(T[1,3]), float(T[2,3]))
        vtk_shape.Update()
        _vtkmesh = vtk_shape.GetOutput()
        # triangle filter
        triangleMesh = vtk.vtkTriangleFilter()  # @UndefinedVariable
        triangleMesh.SetInputData(_vtkmesh)
        triangleMesh.Update()
        triangleMeshOutput = triangleMesh.GetOutput()
        vtkmesh.SetPoints(triangleMeshOutput.GetPoints())
        vtkmesh.SetPolys(triangleMeshOutput.GetPolys()) 
        return vtkmesh
    def slice(self):
        return VTKContours(self)
    def render(self, renderer):
        """
        Render the mesh
         
        :param renderer: a renderer into which an actor derived from this mesh will be added
        :type renderer: `vtk.vtkRenderer`
        """
        assert isinstance(renderer, vtk.vtkRenderer)  # @UndefinedVariable
         
        self.mapper = vtk.vtkOpenGLPolyDataMapper()  # @UndefinedVariable
        # normals
        normals = vtk.vtkPolyDataNormals()  # @UndefinedVariable
        normals.SetInputData(self)
        normals.ComputePointNormalsOn()
        normals.SplittingOff()
        normals.AutoOrientNormalsOn()
        normals.ConsistencyOn()
        normals.Update()
        self.mapper.SetInputData(normals.GetOutput())
        self.actor = vtk.vtkOpenGLActor() # @UndefinedVariable
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
        if self.primary_descriptor == "meshList":
            for mesh in self._sff_segment.meshes:
                self._vtk_meshes.append(
                    VTKMesh.from_mesh(mesh, self.colour, self._vtk_args)
                    )
        elif self.primary_descriptor == "contourList":
            contours = self._sff_segment.contours
            if contours:
                self._vtk_meshes.append(
                    VTKMesh.from_contours(contours, self.colour, self._vtk_args)
                    )
        elif self.primary_descriptor == "shapePrimitiveList":
            for shape in self._sff_segment.shapes:
                transform = self._transforms[shape.transformId]
                self._vtk_meshes.append(
                    VTKMesh.from_shape(shape, self.colour, self._vtk_args, transform)
                    )
        elif self.primary_descriptor == "threeDVolume":
            mask = self._sff_segment.mask
            self._vtk_meshes.append(
                VTKMesh.from_volume(mask, self.colour, self._vtk_args)
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
            self.__setattr__('{}_plane'.format(d), vtk.vtkPlane())  # @UndefinedVariable
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
        filledHoles = vtk.vtkFillHolesFilter()  # @UndefinedVariable
        filledHoles.SetInputData(self._mesh)
        filledHoles.SetHoleSize(1000)
        filledHoles.Update()
        filledHolesOutput = filledHoles.GetOutput()
        # first generate normals
        self._normals = vtk.vtkPolyDataNormals()  # @UndefinedVariable
        self._normals.SetInputData(filledHolesOutput)
        self._normals.Update()
        cut_meshes = dict()
        for d in self._axes:
            """
            :TODO: parallelise
            """
            # the configure the cutter
            self._cutter = vtk.vtkCutter()  # @UndefinedVariable
            self._cutter.SetInputConnection(self._normals.GetOutputPort())
            self._cutter.SetCutFunction(self.__getattribute__('%s_plane' % d))
            # determine the number of z-contours to generate from the bounds
            d_max = math.floor(self.__getattribute__('%smax' % d.upper()))
            d_count = int(2*d_max + 1)
            # generate the cuts on the mesh
            self._cutter.GenerateValues(d_count, -d_max, d_max)
            # pass through vtk.vtkStripper to ensure we have triangel strips
            self._stripper = vtk.vtkStripper()  # @UndefinedVariable
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
        assert isinstance(renderer, vtk.vtkRenderer)  # @UndefinedVariable
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
            self.mapper = vtk.vtkOpenGLPolyDataMapper()  # @UndefinedVariable
            self.mapper.SetInputData(getattr(self, '{}_vtkcontours'.format(d)))
            self.actor = vtk.vtkOpenGLActor() # @UndefinedVariable
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
        colour = list(self._sff_segment.colour.rgba.value)
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
        self._sff_seg = sff_seg # the EMDB-SFF segmentation
        if not args.primary_descriptor:
            args.primary_descriptor = self._sff_seg.primaryDescriptor
        self._vtk_args = args
        self.configs = configs
        self._header = VTKHeader(self._sff_seg, self._vtk_args)
        self._segments = list()
        self._sliced_segments = list()
        if self._vtk_args.primary_descriptor == "threeDVolume":
            # if this is segmentation has threeDVolumes...
            # we need to build a mask for each segment
            volume = self._sff_seg.segments[0].volume
            with h5py.File(os.path.join(self._sff_seg.filePath, volume.file)) as f:
                try:
                    # segger HDF5
                    r_ids = f['/region_ids'].value # placed at the top to fail fast if an EMDB map
                    p_ids = f['/parent_ids'].value
                    mask = f['/mask'].value
                    c_ids = f['/region_colors'].value
                    r_p_zip = zip(r_ids, p_ids)
                    r_c_zip = dict(zip(r_ids, c_ids))
                    simplified_mask, segment_ids = simplify_mask(mask, r_ids, r_p_zip, r_c_zip)
                    for segment in self._sff_seg.segments:
                        if segment.id in segment_ids:
                            new_simplified_mask = numpy.ndarray(mask.shape, dtype=int)  # @UndefinedVariable
                            new_simplified_mask = 0
                            new_simplified_mask = (simplified_mask == segment.id) * segment.id
                            segment.mask = new_simplified_mask # set the empty 'mask' attribute
                            self._segments.append(VTKSegment(segment, self._vtk_args))
                except KeyError: # EMDB map HDF5
                    for segment in self._sff_seg.segments:
                        segment.mask = f['/mask'].value
                        self._segments.append(VTKSegment(segment, self._vtk_args))
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
    
    def slice(self):
        for segment in self.segments:
            segment.slice()
    
    def as_roi(self, configs):
        from ..formats.roi import ROISegmentation
        return ROISegmentation.from_vtk(self, configs)
    
    def render(self):
        """Render to display"""
        # define the renderer                
        ren = vtk.vtkOpenGLRenderer()  # @UndefinedVariable
        ren.SetBackground(*self._vtk_args.background_colour)
        # populate the renderer with the meshes
        for segment in self.segments:
            ren = segment.render(ren)
        # render window
        renWin = vtk.vtkRenderWindow()  # @UndefinedVariable
        renWin.AddRenderer(ren)
        if self._vtk_args.full_screen:
            renWin.FullScreenOn()
        # render window interactor
        iren = vtk.vtkRenderWindowInteractor()  # @UndefinedVariable
        iren.SetRenderWindow(renWin)
        iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())  # @UndefinedVariable
        # from vtkCubeAxesActor.h
        #define VTK_FLY_OUTER_EDGES     0
        #define VTK_FLY_CLOSEST_TRIAD   1
        #define VTK_FLY_FURTHEST_TRIAD  2
        #define VTK_FLY_STATIC_TRIAD    3
        #define VTK_FLY_STATIC_EDGES    4
        if self._vtk_args.cube_axes is not None:
            cubeAxesActor = vtk.vtkCubeAxesActor()  # @UndefinedVariable
            cubeAxesActor.SetBounds(ren.ComputeVisiblePropBounds())
            cubeAxesActor.SetCamera(ren.GetActiveCamera())
            cubeAxesActor.SetFlyMode(self._vtk_args.cube_axes)
            cubeAxesActor.SetFlyModeToStaticEdges() # how the cube axes will appear
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
            axesActor = vtk.vtkAxesActor()  # @UndefinedVariable
            axesWidget = vtk.vtkOrientationMarkerWidget()  # @UndefinedVariable
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
        for segment in self.segments:
            segment_data = dict()
            segment_data['id'] = segment.id
            segment_data['colour'] = map(float, segment.colour)
            segment_data['meshes'] = list()
            for j, mesh in enumerate(segment.meshes):
                decimateMesh = vtk.vtkDecimatePro()  # @UndefinedVariable
                decimateMesh.SetInputData(mesh)
                decimateMesh.SetTargetReduction(0.9)
                decimateMesh.PreserveTopologyOn()
                decimateMesh.Update()
                decimateMeshOutput = decimateMesh.GetOutput()
                writer = vtk.vtkXMLPolyDataWriter()  # @UndefinedVariable
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