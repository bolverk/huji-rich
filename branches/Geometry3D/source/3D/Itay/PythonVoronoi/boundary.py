"""
Implementation of boundary related classes - used to copy and reflect points,
and check for hull breaches
"""
__author__ = 'zmbq'
import numpy as np
import numpy.ma as ma
from numpy import linalg

class Subcube:
    """
    A Subcube describes one of the 27 subcubes of a big cube. A subcube is determined by three
    offsets, '+', '-' or ' ' for each axis.
    For example, "   " is the main subcube (at (0,0,0)), "+  " is the subcube to the right of the
    main subcube. "--+" is the subcube to the left, below and in front of the main subcube.
    """
    def __init__(self, offsets):
        self._check_offsets(offsets)
        self._offsets = tuple(offsets)

    def _check_offsets(self, offsets):
        """
        Makes sure offsets is a legal offset string
        :return: Throws an exception if it isn't.
        """
        if len(offsets)!=3:
            raise ValueError("Subcube offsets must be 3 characters long")
        allowed = [' ', '-', '+']
        for offset in offsets:
            if not offset in allowed:
                raise ValueError("Invalid offest '%s'" % offsets)
        return  # All is well

    def __str__(self):
        return self._offsets

    def __getitem__(self, index):
        return self._offsets[index]

    @property
    def num_offseted(self):
        """
        :return: The number of 'offseted' axes. That is, axes that don't have a ' ' as their offset
        """
        return sum([1 if offset!=' ' else 0 for offset in self._offsets])

class Boundary:
    def __init__(self, ftr, bll):
        """
        Constructs a Boundary
        :param: tfr The Front Top Right corner of the boundary
        :param: bll The Back Lower Left corner of the boundary
        """
        if len(ftr)!=3 or len(bll)!=3:
            raise ValueError("3D corners are expected")
        for i in range(3):
            if ftr[i] <= bll[i]:
                raise ValueError("tfr must be in front of, above and to the right of bll")
        self._ftr = ftr
        self._bll = bll

    @property
    def ftr(self):
        return self._ftr

    @property
    def bll(self):
        return self._bll

    def __contains__(self, item):
        return self.is_inside(item)

    def is_inside(self, pt):
        """
        Checks if a point is inside the boundary
        """
        if len(pt)!=3:
            raise ValueError("3D points expected")
        for axis in range(3):
            if not self._bll[axis] <= pt[axis] <= self._ftr[axis]:
                return False
        return True

def vector_to_subcube(point, boundary, subcube):
    """
    Returns the shortest vector from point to the subcube
    :param point: The point
    :param boundary: The boundary
    :param subcube: The subcube
    :return: The shortest vector. The vector's length if the distance to the subcube
    """

    def subcube_mask(boundary, subcube):
        """
        Returns a masked array with the coordinates on the boundary depending on the subcube.
        For each + offset, the ftr coordinate is taken, for each -, the bll coordinate is
        taken. For ' ', the value is masked
        """
        pt = np.zeros(3, dtype=np.float64)
        mask = np.ones(3, dtype=np.bool)
        for axis in range(3):
            if subcube[axis]=='-':
                pt[axis] = boundary.bll[axis]
                mask[axis] = False
            elif subcube[axis]=='+':
                pt[axis] = boundary.ftr[axis]
                mask[axis] = False

        return ma.masked_array(pt, mask, fill_value=0)

    masked = subcube_mask(boundary, subcube)
    return ma.filled(masked - point)

def distance_from_subcube(point, boundary, subcube):
    vec = vector_to_subcube(point, boundary, subcube)
    return linalg.norm(vec)

def reflect_points(pts, boundary, offsets):
    """
    Reflect points through the boundary, using only the supplied subcubes
    """
    reflected = np.zeros([len(pts) * len(offsets), 3], dtype=np.float64)
    num = 0
    subcubes = [Subcube(offset) for offset in offsets]
    for pt in pts:
        for subcube in subcubes:
            vec = vector_to_subcube(pt, boundary, subcube)
            reflected[num] = pt + 2*vec
            num += 1

    return reflected