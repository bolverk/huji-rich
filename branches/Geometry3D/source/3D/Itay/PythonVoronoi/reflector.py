from TetGenReader import read_tetgen_points
from boundary import Boundary, reflect_point, Subcube, distance_from_subcube
import numpy as np

__author__ = 'zmbq'

import os.path

def _subcubes_reflector(pts, boundary, subcubes):
    reflected = np.zeros((len(pts) * len(subcubes), 3), dtype=np.float64)

    num = 0
    for i in range(len(pts)):
        pt = pts[i]
        for subcube in subcubes:
            reflected[num] = reflect_point(pt, boundary, subcube)
            num += 1

    return reflected

_brute_force_offsets = ('-  ', ' - ', '  -', '  +', ' + ', '+  ')
_brute_force_subcubes = [Subcube(o) for o in _brute_force_offsets]
def brute_force_reflector(pts, boundary):
    return  _subcubes_reflector(pts, boundary, _brute_force_subcubes)

_all_offsets = ("---", "-- ", "--+", "- -", "-  ", "- +", "-+-", "-+ ", "-++",
                " --", " - ", " -+", "  -",        "  +", " +-", " + ", " ++",
                "+--", "+- ", "+-+", "+ -", "+  ", "+ +", "++-", "++ ", "+++",)
_all_subcubes = [Subcube(o) for o in _all_offsets]
def full_brute_force_reflector(pts, boundary):
    return _subcubes_reflector(pts, boundary, _all_subcubes)

class CloseToBoundaryReflector:
    def __init__(self, points, tetrahedra, centers, radii, tetrahedron_neighbors, edge_neighbors, vertex_neighbors):
        self._points = points
        self._tetrahedra = tetrahedra
        self._centers = centers
        self._radii = radii
        self._tetrahedron_neighbors = tetrahedron_neighbors
        self._edge_neighbors = edge_neighbors
        self._vertex_neighbors = vertex_neighbors

        self._outer_tetrahedra = self._find_outer()
        self._edge_tetrahedra = self._find_edge()

    def _is_big_tetrahedron(self, pt_num):
        return pt_num >= len(self._points) - 4

    def _find_outer(self):
        """
        Find the outer tetrahedra - tetrahedra that touch the big tetrahedron
        """
        outer = []
        for i in range(len(self._tetrahedra)):
            for pt in self._tetrahedra[i]:
                if self._is_big_tetrahedron(pt):
                    outer.append(i)
                    break

        return outer

    def _find_edge(self):
        """
        Find the edge tetrahedra - tetrahedra that share a vertex with an outer tetrahedra
        """
        edge = []
        for i in range(len(self._tetrahedra)):
            if i in self._outer_tetrahedra:
                continue
            tetrahedron = self._tetrahedra[i]
            is_edge = False
            for pt in tetrahedron:
                for neighbor in self._vertex_neighbors[pt]:
                    if neighbor in self._outer_tetrahedra:
                        edge.append(i)
                        is_edge = True
                        break
                if is_edge:
                    break
        return edge

    def _hull_breaches(self, boundary):
        """
        Returns the hull breaches of all the points: {pt->[subcube offsets]}
        """
        breaches = {}
        candidates = self._edge_tetrahedra

        while candidates:
            next_candidates = set()
            for tetrahedron in candidates:
                for pt in self._tetrahedra[tetrahedron]:
                    if not pt in breaches:
                        breaches[pt] = []
                    pt_neighbors = self._vertex_neighbors[pt]
                    for touching_tetrahedron in pt_neighbors:
                        if touching_tetrahedron in self._outer_tetrahedra:
                            continue
                        center = self._centers[touching_tetrahedron]
                        radius = self._radii[touching_tetrahedron]

                        breaching = False
                        for subcube in _brute_force_subcubes:
                            distance = distance_from_subcube(center, boundary, subcube)
                            if distance < radius: # Breaching tetrahedron!
                                if not subcube.offsets in breaches[pt]: # A new breach!!
                                    breaches[pt].append(subcube.offsets)
                                    breaching = True

                        if breaching:
                            next_candidates.add(touching_tetrahedron)
                            next_candidates.union(set(pt_neighbors))
            if -1 in next_candidates:
                next_candidates.remove(-1)
            candidates = next_candidates
        return breaches

    def __call__(self, pts, boundary):
        breaches = self._hull_breaches(boundary)
        num_breaches = sum([len(breaches[breach]) for breach in breaches])
        reflected = np.zeros((num_breaches, 3), dtype=np.float64)

        num = 0
        for ptidx in breaches.keys():
            pt = self._points[ptidx]
            for offset in breaches[ptidx]:
                subcube = Subcube(offset)
                reflected[num] = reflect_point(pt, boundary, subcube)
                num += 1

        return reflected