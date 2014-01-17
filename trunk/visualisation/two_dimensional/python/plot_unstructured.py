#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys
import math

class Vector2D:
    """
    A mathematical vector in 2d
    """

    def __init__(self, x, y):
        """
        Class constructor
        Input:
        x - x Component
        y - y Component
        """

        self._x = x
        self._y = y
        return

    def SetX(x):
        """
        Sets the x component of a vector
        Input:
        x - x Component
        """

        self._x = x
        return

    def SetY(y):
        """
        Sets the y component of a vector
        Input:
        y - y Component
        """

        self._y = y
        return
    
    def GetX(self):
        """
        Returns the x component of a vector
        """

        return self._x

    def GetY(self):
        """
        Returns the y component of a vector
        """

        return self._y

    def __add__(self, v):
        """
        Adds the components of one vector to another
        """

        return Vector2D(self.GetX()+v.GetX(), \
                            self.GetY()+v.GetY())

    def __sub__(self, v):
        """
        Subtracting the components of one vector from another
        """

        return Vector2D(self.GetX()-v.GetX(), \
                            self.GetY()-v.GetY())

    def __mul__(self, c):
        """
        Multiplying each vector term by a scalar
        """

        return Vector2D(self.GetX()*c, self.GetY()*c)

    def __div__(self, c):
        """
        Dividing each vector term by a scalar
        """

        return Vector2D(self.GetX()/c, self.GetY()/c)

def crossz(v):
    
    return Vector2D(v.GetY(), -v.GetX())

class Edge:
    """
    Interface between two cells
    """

    def __init__(self, vertices, neighbors):
        """
        Class constructor
        Input:
        vertices - List of vertices
        neighbors - List of the indices of neighbor cells
        """
        self._vertices = vertices
        self._neighbors = neighbors
        return

    def GetVertices(self):
        """
        Returns a list of vertices
        """

        return self._vertices

    def GetNeighbors(self):
        """
        Returns a list of indices of neibor cells
        """

        return self._neighbors

def ScalarProd(v1, v2):
    
    return v1.GetX()*v2.GetX()+ \
        v1.GetY()*v2.GetY()

def Modulus(v):

    return ScalarProd(v,v)**0.5

def SubtendedAngle(v1, v2):

    aux = ScalarProd(v1,v2)/(Modulus(v1)*Modulus(v2))
    if(aux>1):
        aux = 1
    if(aux<-1):
        aux = -1
    return math.acos(aux)

def CalcCentroid(edge):
    
    vertices = edge.GetVertices()
    return (vertices[0]+vertices[1])*0.5

def CalcCenter(edges):

    center = Vector2D(0,0)
    for i in edges:
        center = center + CalcCentroid(i)/len(edges)
    return center

def  sgn(x):

    if(x>0):
        return 1
    if(x<0):
        return -1
    else:
        return 0

def DirAngle(v, vref):

    return SubtendedAngle(v, vref)*sgn(ScalarProd(v,crossz(vref)))

def IsRightTriangle(p1, p2, p3):

    return ScalarProd(p3-p1,crossz(p2-p1))>=0

class Polygon:
    """
    Convex polygon
    """

    def __init__(self, edges):
        """
        Class constructor
        Input:
        edges - Unordered list of vertices
        """
        
        self._center = CalcCenter(edges)
        self._vertices = []
        for i in edges:
            edge_vertices = i.GetVertices()
            if IsRightTriangle(self._center, edge_vertices[0], edge_vertices[1]):
                self._vertices.append(edge_vertices[0])
            else:
                self._vertices.append(edge_vertices[1])
        for i in range(len(self._vertices)):
            for j in range(i+1,len(self._vertices)):
                for k in range(j+1,len(self._vertices)):
                    if ~IsRightTriangle(self._vertices[i], \
                                            self._vertices[j], \
                                            self._vertices[k]):
                        self._vertices[j], self._vertices[k] = self._vertices[k], self._vertices[j]
    
    def GetVertices(self):
        """
        Returns a list of vertices, in a clockwise order
        """

        return self._vertices

    def GetCenter(self):
        """
        Returns the center of the polygon
        """

        return self._center



# Code starts here

rawd = numpy.loadtxt(sys.argv[1])
z = numpy.loadtxt(sys.argv[2])
zeps = 1e-15
zn = (z-z.min())/(z.max()-z.min()+zeps)

# Create edges
edges = []
for i in rawd:
    v1 = Vector2D(i[0],i[1])
    v2 = Vector2D(i[3],i[4])
    i1 = int(i[2])
    i2 = int(i[5])
    e0 = Edge([v1,v2],[i1,i2])
    edges.append(e0)

# Sort edges by polygon
upe = [] # unordered polygon edges
for i in range(len(z)):
    upe.append([])
for i in edges:
    for j in i.GetNeighbors():
        if j>-1:
            upe[j].append(i)

counter = 0
for i,j in zip(upe, zn):

    aux = Polygon(i)
    vert_x_list = [k.GetX() for k in aux.GetVertices()]
    vert_y_list = [k.GetY() for k in aux.GetVertices()]
    pylab.fill(vert_x_list, vert_y_list,
               fc=matplotlib.cm.jet(j))

pylab.axis('equal')

ax, _ = matplotlib.colorbar.make_axes(pylab.gca())
matplotlib.colorbar.ColorbarBase(ax,norm=matplotlib.colors.Normalize(vmin=min(z),vmax=max(z)))

pylab.show()
