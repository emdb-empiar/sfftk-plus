# -*- coding: utf-8 -*-
from __future__ import division
"""
roi_primitives.py
========================

Module containing wrapper classes to work with ROI objects


Requirements
==================================
- the path to Python library that ships with OMERO.server


References
==========
* Model Overview: http://www.openmicroscopy.org/site/support/ome-model/developers/model-overview.html
* OMERO Python Examples: https://www.openmicroscopy.org/site/support/omero5.2/developers/Python.html
* OMERO Python Bindings: http://downloads.openmicroscopy.org/omero/5.2.1/api/python/index.html
* OMERO Python Model Objects: http://downloads.openmicroscopy.org/omero/5.2.1/api/python/omero/omero.model.html
* ROI Model: http://www.openmicroscopy.org/site/support/ome-model/developers/roi.html
* OME Schema: http://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2015-01/ome.html
* BlitzGateway Docs: https://www.openmicroscopy.org/site/support/omero5.2/developers/PythonBlitzGateway.html


TODO:
- documentation
- 


Version history:
0.0.1, 2016-01-20, First incarnation
0.0.2, 2016-02-22, Moved and renamed
0.0.3, 2016-06-02, Removed ROI class

"""

__author__ 	= 'Paul K. Korir, PhD'
__email__ 	= 'pkorir@ebi.ac.uk, paul.korir@gmail.com'
__date__ 	= '2016-01-20'

import sys
import omero.model
from omero.rtypes import rint, rdouble, rstring
from omero.model.enums import UnitsLength


from omero_version import omero_version

class Shape(omero.model.Shape):
	def __init__(self, theZ=None, theT=None, theC=None, *args, **kwargs):
		super(Shape, self).__init__(*args, **kwargs)
		if theZ is not None:
			self.setTheZ(theZ)
		if theT is not None:
			self.setTheT(theT)
		if theC is not None:
			self.setTheC(theC)

	def getFillColor(self, *args, **kwargs):
		return super(Shape, self).getFillColor(*args, **kwargs).getValue()

	def setFillColor(self, R, G, B, A=1.0, normalised=True, *args, **kwargs):
		if normalised:
			R = int(255*R)
			G = int(255*G)
			B = int(255*B)
			if A == 0:
				print >> sys.stderr, 'Setting alpha to 1% because 0% does not work'
				A = 1/255
			A = int(255*A)
# 			else:
# 				A = 0
# 		rgba_int = (A << 24) + (R << 16) + (G << 8) + B
		rgba_int = (R << 24) + (G << 16) + (B << 8) + A
		if rgba_int > 2147483647:
			rgba_int = -2147483648 + rgba_int % 2147483648
		theFillColor = rint(rgba_int)
		super(Shape, self).setFillColor(theFillColor, *args, **kwargs)

	def getFontFamily(self, *args, **kwargs):
		return super(Shape, self).getFontFamily(*args, **kwargs).getValue()

	def setFontFamily(self, theFontFamily, *args, **kwargs):
		"""Set the font family.

		theFontFamily
			'serif', 'san-serif', 'cursive', 'fantasy', 'monospace'
		"""
		fontFamilies = [
			'serif', 
			'sans-serif', 
			'cursive',
			'fantasy',
			'monospace'
			]
		try:
			assert theFontFamily in fontFamilies
			super(Shape, self).setFontFamily(rstring(theFontFamily), *args, **kwargs)
		except:
			raise ValueError('unknown font family: %s' % theFontFamily)

	def getFontSize(self, *args, **kwargs):
		return super(Shape, self).getFontSize(*args, **kwargs).getValue()

	def setFontSize(self, theFontSize, unit, *args, **kwargs):
		theFontSize = omero.model.LengthI(theFontSize, UnitsLength.__dict__[unit.upper()])
		super(Shape, self).setFontSize(theFontSize, *args, **kwargs)

	def getFontStyle(self, *args, **kwargs):
		return super(Shape, self).getFontStyle(*args, **kwargs).getValue()

	def setFontStyle(self, fontStyle, *args, **kwargs):
		fontStyle = rstring(fontStyle)
		super(Shape, self).setFontStyle(fontStyle, *args, **kwargs)

	def getStrokeColor(self, *args, **kwargs):
		return super(Shape, self).getStrokeColor(*args, **kwargs).getValue()

	def setStrokeColor(self, R, G, B, A=1.0, normalised=True, *args, **kwargs):
		if normalised:
			R = int(255*R)
			G = int(255*G)
			B = int(255*B)
			if A == 0:
				print >> sys.stderr, 'Setting alpha to 1% because 0% does not work'
				A = 1/255
			A = int(255*A)
# 		rgba_int = (A << 24) + (R << 16) + (G << 8) + B
		rgba_int = (R << 24) + (G << 16) + (B << 8) + A
		if rgba_int > 2147483647:
			rgba_int = -2147483648 + rgba_int % 2147483648
		theStrokeColor = rint(rgba_int)
		super(Shape, self).setStrokeColor(theStrokeColor, *args, **kwargs)

	def getStrokeWidth(self, *args, **kwargs):
		return super(Shape, self).getStrokeWidth(*args, **kwargs).getValue()

	def setStrokeWidth(self, theStrokeWidth, unit, *args, **kwargs):
		theStrokeWidth = omero.model.LengthI(theStrokeWidth, UnitsLength.__dict__[unit.upper()])
		super(Shape, self).setStrokeWidth(theStrokeWidth, *args, **kwargs)

	def getTheC(self, *args, **kwargs):
		return super(Shape, self).getTheC(*args, **kwargs).getValue()

	def setTheC(self, theTheC, *args, **kwargs):
		theTheC = rint(theTheC)
		super(Shape, self).setTheC(theTheC, *args, **kwargs)

	def getTheT(self, *args, **kwargs):
		return super(Shape, self).getTheT(*args, **kwargs).getValue()

	def setTheT(self, theTheT, *args, **kwargs):
		theTheT = rint(theTheT)
		super(Shape, self).setTheT(theTheT, *args, **kwargs)

	def getTheZ(self, *args, **kwargs):
		return super(Shape, self).getTheZ(*args, **kwargs).getValue()

	def setTheZ(self, theTheZ, *args, **kwargs):
		theTheZ = rint(theTheZ)
		super(Shape, self).setTheZ(theTheZ, *args, **kwargs)

	def getTextValue(self, *args, **kwargs):
		return super(Shape, self).getTextValue(*args, **kwargs).getValue()

	def setTextValue(self, textValue, *args, **kwargs):
		textValue = rstring(textValue)
		super(Shape, self).setTextValue(textValue, *args, **kwargs)

	# def getTransform(self, *args, **kwargs):
	# 	return super(Shape, self).getTransform(*args, **kwargs)

	# def setTransform(self, A00, A01, A02, A10, A11, A12, *args, **kwargs):
	# 	theTransform = [
	# 		rfloat(A00), rfloat(A01), rfloat(A02),
	# 		rfloat(A10), rfloat(A11), rfloat(A12),
	# 		]
	# 	super(Shape, self).setTransform(theTransform, *args, **kwargs)


class Line(Shape, omero.model.LineI):
	def __init__(self, X1=None, X2=None, Y1=None, Y2=None, *args, **kwargs):
		super(Line, self).__init__(*args, **kwargs)
		# all() is True if the list is empty
		# the list is empty is X1, X2, Y1, Y2 are numbers that can be coerced to rdouble
		if all(filter(lambda x: isinstance(x, int) or isinstance(x, float), [X1, X2, Y1, Y2])):
			self.setX1(X1)
			self.setX2(X2)
			self.setY1(Y1)
			self.setY2(Y2)

	def getX1(self, *args, **kwargs):
		return super(Line, self).getX1(*args, **kwargs).getValue()

	def setX1(self, X1, *args, **kwargs):
		X1 = rdouble(X1)
		super(Line, self).setX1(X1, *args, **kwargs)

	def getX2(self, *args, **kwargs):
		return super(Line, self).getX2(*args, **kwargs).getValue()

	def setX2(self, X2, *args, **kwargs):
		X2 = rdouble(X2)
		super(Line, self).setX2(X2, *args, **kwargs)

	def getY1(self, *args, **kwargs):
		return super(Line, self).getY1(*args, **kwargs).getValue()

	def setY1(self, Y1, *args, **kwargs):
		Y1 = rdouble(Y1)
		super(Line, self).setY1(Y1, *args, **kwargs)

	def getY2(self, *args, **kwargs):
		return super(Line, self).getY2(*args, **kwargs).getValue()

	def setY2(self, Y2, *args, **kwargs):
		Y2 = rdouble(Y2)
		super(Line, self).setY2(Y2, *args, **kwargs)


class Rect(object):
	def __init__(self, *args, **kwargs):
		super(Rect, self).__init__(*args, **kwargs)

	def getX(self, *args, **kwargs):
		return super(Rect, self).getX(*args, **kwargs).getValue()

	def setX(self, X, *args, **kwargs):
		X = rdouble(X)
		super(Rect, self).setX(X, *args, **kwargs)

	def getY(self, *args, **kwargs):
		return super(Rect, self).getY(*args, **kwargs).getValue()

	def setY(self, Y, *args, **kwargs):
		Y = rdouble(Y)
		super(Rect, self).setY(Y, *args, **kwargs)

	def getWidth(self, *args, **kwargs):
		return super(Rect, self).getWidth(*args, **kwargs).getValue()

	def setWidth(self, Width, *args, **kwargs):
		Width = rdouble(Width)
		super(Rect, self).setWidth(Width, *args, **kwargs)

	def getHeight(self, *args, **kwargs):
		return super(Rect, self).getHeight(*args, **kwargs).getValue()

	def setHeight(self, Height, *args, **kwargs):
		Height = rdouble(Height)
		super(Rect, self).setHeight(Height, *args, **kwargs)


class Rectangle(Shape, Rect, omero.model.RectangleI):
	def __init__(self, X=None, Y=None, Width=None, Height=None, *args, **kwargs):
		super(Rectangle, self).__init__(*args, **kwargs)
		# all() is True if the list is empty
		# the list is empty is X, Y, Width, Height are numbers that can be coerced to rdouble
		if all(filter(lambda x: isinstance(x, int) or isinstance(x, float), [X, Y, Width, Height])):
			self.setX(X)
			self.setY(Y)
			self.setWidth(Width)
			self.setHeight(Height)


class Mask(Shape, Rect, omero.model.MaskI):
	def __init__(self, X=None, Y=None, Width=None, Height=None, *args, **kwargs):
		super(Mask, self).__init__(*args, **kwargs)
		# all() is True if the list is empty
		# the list is empty is X, Y, Width, Height are numbers that can be coerced to rdouble
		if all(filter(lambda x: isinstance(x, int) or isinstance(x, float), [X, Y, Width, Height])):
			self.setX(X)
			self.setY(Y)
			self.setWidth(Width)
			self.setHeight(Height)


class Ellipse(Shape, omero.model.EllipseI):
	def __init__(self, X=None, Y=None, RadiusX=None, RadiusY=None, *args, **kwargs):
		super(Ellipse, self).__init__(*args, **kwargs)
		# all() is True if the list is empty
		# the list is empty is X, Y, Rx, Ry are numbers that can be coerced to rdouble
		if all(filter(lambda x: isinstance(x, int) or isinstance(x, float), [X, Y, RadiusX, RadiusY])):
			self.setX(X)
			self.setY(Y)
			self.setRadiusX(RadiusX)
			self.setRadiusY(RadiusY)

	def getX(self, *args, **kwargs):
		return super(Ellipse, self).getX(*args, **kwargs).getValue()

	def setX(self, Cx, *args, **kwargs):
		super(Ellipse, self).setX(rdouble(Cx), *args, **kwargs)

	def getY(self, *args, **kwargs):
		return super(Ellipse, self).getY(*args, **kwargs).getValue()

	def setY(self, Cy, *args, **kwargs):
		super(Ellipse, self).setY(rdouble(Cy), *args, **kwargs)

	def getRadiusX(self, *args, **kwargs):
		return super(Ellipse, self).getRadiusX(*args, **kwargs).getValue()

	def setRadiusX(self, Rx, *args, **kwargs):
		super(Ellipse, self).setRadiusX(rdouble(Rx), *args, **kwargs)

	def getRadiusY(self, *args, **kwargs):
		return super(Ellipse, self).getRadiusY(*args, **kwargs).getValue()

	def setRadiusY(self, Ry, *args, **kwargs):
		super(Ellipse, self).setRadiusY(rdouble(Ry), *args, **kwargs)


class Point(Shape, omero.model.PointI):
	def __init__(self, X=None, Y=None, *args, **kwargs):
		super(Point, self).__init__(*args, **kwargs)
		# all() is True if the list is empty
		# the list is empty is Cx, Cy are numbers that can be coerced to rdouble
		if all(filter(lambda x: isinstance(x, int) or isinstance(x, float), [X, Y])):
			self.setX(X)
			self.setY(Y)

	def getX(self, *args, **kwargs):
		return super(Point, self).getX(*args, **kwargs).getValue()

	def setX(self, X, *args, **kwargs):
		super(Point, self).setX(rdouble(X), *args, **kwargs)

	def getY(self, *args, **kwargs):
		return super(Point, self).getY(*args, **kwargs).getValue()

	def setY(self, Y, *args, **kwargs):
		super(Point, self).setY(rdouble(Y), *args, **kwargs)


class Poly(Shape):
	def __init__(self, *args, **kwargs):
		super(Poly, self).__init__(*args, **kwargs)

	def getPoints(self, *args, **kwargs):
		points = super(Poly, self).getPoints(*args, **kwargs)
		if points is not None:
			import re
			m = re.match(r'points\[(?P<points0>[0-9., \-]+)\] points1\[(?P<points1>[0-9., \-]+)\] points2\[(?P<points2>[0-9., \-]+)\] mask\[0,0,0\]', points.getValue())
			points0 = map(lambda p: tuple(map(float, p.split(','))), m.group('points0').split(', '))
			return points0
		else:
			return []

	def setPoints(self, points, offsetXFrom=None, offsetYFrom=None, swapXY=False, *args, **kwargs):
		"""Construct a Shape from points.

		points
			vector of 2-tuples (x,y)
		offsetXFrom
			positive integer. Use offsetXFrom-x instead of x
		offsetYFrom
			positive integer. Use offsetYFrom-y instead of y
		swapXY
			use (y,x) instead of (x,y)
		"""
		if offsetXFrom is not None:
			assert offsetXFrom > 0
			assert isinstance(offsetXFrom, int)
		if offsetYFrom is not None:
			assert offsetYFrom > 0
			assert isinstance(offsetYFrom, int)
			
		if swapXY:
			if offsetXFrom is not None:
				if offsetYFrom is not None:
					points = [(offsetYFrom-y,offsetXFrom-x) for x,y in points]
				else:
					points = [(y,offsetXFrom-x) for x,y in points]
			else:
				if offsetYFrom is not None:
					points = [(offsetYFrom-y,x) for x,y in points]
				else:
					points = [(y,x) for x,y in points]
		else:
			if offsetXFrom is not None:
				if offsetYFrom is not None:
					points = [(offsetXFrom-x,offsetYFrom-y) for x,y in points]
				else:
					points = [(offsetXFrom-x,y) for x,y in points]
			else:
				if offsetYFrom is not None:
					points = [(x,offsetYFrom-y) for x,y in points]
				else:
					points = [(x,y) for x,y in points]

		p = ", ".join([str(a) + "," + str(b) for a,b in points])
		m = ",".join(map(str, [0]*3))
		points_string = rstring("points[%s] points1[%s] points2[%s] mask[%s]" % (p, p, p, m))

		# finally, call the super
		super(Poly, self).setPoints(points_string, *args, **kwargs)


class Polyline(Poly, omero.model.PolylineI):
	def __init__(self, *args, **kwargs):
		super(Polyline, self).__init__(*args, **kwargs)


class Polygon(Poly, omero.model.PolygonI):
	def __init__(self, *args, **kwargs):
		super(Polygon, self).__init__(*args, **kwargs)


class Label(Shape, omero.model.LabelI):
	def __init__(self, X=None, Y=None, label=None, *args, **kwargs):
		super(Label, self).__init__(*args, **kwargs)
		# all() is True if the list is empty
		# the list is empty is X, Y are numbers that can be coerced to rdouble
		if all(filter(lambda x: isinstance(x, int) or isinstance(x, float), [X, Y])):
			self.setX(X)
			self.setY(Y)
		if label is not None:
			self.setTextValue(label)

	def getX(self, *args, **kwargs):
		return super(Label, self).getX(*args, **kwargs).getValue()

	def setX(self, X, *args, **kwargs):
		super(Label, self).setX(rdouble(X), *args, **kwargs)

	def getY(self, *args, **kwargs):
		return super(Label, self).getY(*args, **kwargs).getValue()

	def setY(self, Y, *args, **kwargs):
		super(Label, self).setY(rdouble(Y), *args, **kwargs)
