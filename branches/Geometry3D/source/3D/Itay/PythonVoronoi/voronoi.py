__author__ = 'zmbq'

import numpy as np
import numpy.linalg as linalg
import itertools
from TetGenReader import read_tetgen

orig_points = None
all_points = None
tetrahedra = None
centers = None
neighbors = None
boundary = ((-200, -200, -200), (200, 200, 200))

def svec(vec):
    return "(%.5f, %.5f, %.5f)" % (vec[0], vec[1], vec[2])

names = {(0,0,0,): 'Z0'}
name_to_vec = {'Z0': (0,0,0,)}

def get_name(vector, prefix='V'):
    vector = tuple(vector)
    for vec in names.keys():
        dist = linalg.norm(vector - np.array(vec)) ** 2
        if dist < 1e-5:
            return names[vec]
    name = prefix + str(len(names))
    names[vector] = name
    name_to_vec[name] = vector
    return name

def is_original_point(vector_idx):
    return vector_idx < len(orig_points)

def tetrahedron_coords(tetrahedron, points):
    """ Converts a tetrahedron from indices to 3D points
    :param tetrahedron: Tetrahedron in indices to points
    :param points: The entire TetGen point array
    :return: An array of the 4 coordinates
    """
    coords = [points[t] for t in tetrahedron]
    return coords

def tetrahedron_center(tetrahedron):
    """
    Calculates the tetrahedron center
    :param tetrahedron: The tetrahedron in 3D coordinates
    :return: The center in 3D coordinates
    Taken from here: http://mathworld.wolfram.com/Circumsphere.html
    """
    v1 = tetrahedron[0]
    v2 = tetrahedron[1]
    v3 = tetrahedron[2]
    v4 = tetrahedron[3]
    m_a = np.array([[v1[0], v1[1], v1[2], 1],
                   [v2[0], v2[1], v2[2], 1],
                   [v3[0], v3[1], v3[2], 1],
                   [v4[0], v4[1], v4[2], 1]])
    a = linalg.det(m_a)

    m_Dx = np.array([[linalg.norm(v1)**2, v1[1], v1[2], 1],
                    [linalg.norm(v2)**2, v2[1], v2[2], 1],
                    [linalg.norm(v3)**2, v3[1], v3[2], 1],
                    [linalg.norm(v4)**2, v4[1], v4[2], 1]])
    Dx = linalg.det(m_Dx)

    m_Dy = np.array([[linalg.norm(v1)**2, v1[0], v1[2], 1],
                    [linalg.norm(v2)**2, v2[0], v2[2], 1],
                    [linalg.norm(v3)**2, v3[0], v3[2], 1],
                    [linalg.norm(v4)**2, v4[0], v4[2], 1]])
    Dy = linalg.det(m_Dy)

    m_Dz = np.array([[linalg.norm(v1)**2, v1[0], v1[1], 1],
                    [linalg.norm(v2)**2, v2[0], v2[1], 1],
                    [linalg.norm(v3)**2, v3[0], v3[1], 1],
                    [linalg.norm(v4)**2, v4[0], v4[1], 1]])
    Dz = linalg.det(m_Dz)

    return np.array([Dx, -Dy, Dz] / (2*a))

def save_point_names():
    for point in orig_points:
        get_name(point, 'C')

def calc_centers():
    tetrahedra_coords = [tetrahedron_coords(t, all_points) for t in tetrahedra]
    centers = [tetrahedron_center(t) for t in tetrahedra_coords]
    return centers

class Face:
    def __init__(self, vertices):
        self._vertices = frozenset(vertices)
        self._neighbors = []

    @property
    def vertices(self):
        return self._vertices

    @property
    def neighbors(self):
        return self._neighbors

    def add_neighbor(self, cell):
        if cell in self._neighbors:
            return
        if len(self._neighbors)==2:
            raise ValueError("Can't add more than two neighbors to a face")
        self._neighbors.append(cell)

    def same_vertices(self, vertices):
        vset = frozenset(vertices)
        return self._vertices == vset
        #if len(vset)!=len(self._vertices):
        #    return False
        #for vertex in vset:
        #    if not vertex in self._vertices:
        #        return False
        #return True

edges = {}  # (vec1, vec2) -> [list of tetrahedra touching edge]
def normalize_edge(edge):
    if edge[0] < edge[1]:
        return edge[1], edge[0]
    return edge

def build_edges():
    edges = {}
    num = 0
    for tetrahedron in tetrahedra:
        for edge in itertools.combinations(tetrahedron, 2):
            n = normalize_edge(edge)
            if not n in edges:
                edges[n] = []
            edges[n].append(num)
        num += 1
    return edges

faces = []
def save_face(vertices):
    for i in range(len(faces)):
        if faces[i].same_vertices(vertices):
            return i
    face = Face(vertices)
    faces.append(face)
    return len(faces) - 1

def create_face(edge):
    edge_neighbors = edges[edge]
    if not edge_neighbors:
        return None
    tetrahedron_centers = [centers[n] for n in edge_neighbors]
    vertex_names = [get_name(c) for c in tetrahedron_centers]
    vertex_set = frozenset(vertex_names)
    if len(vertex_set) < 3:   # Degenerate 0-area face, ignore it
        return None
#    name_set = set(vertex_names)
#    if(len(name_set) < len(vertex_names)):
#        return None  # Repeating vertices - this is a degenerate face

    return save_face(vertex_set)

cells = []

def voronoi():
    global cells
    for tetrahedron in tetrahedra:
        for edge in itertools.combinations(tetrahedron, 2):
            if not is_original_point(edge[0]) and not is_original_point(edge[1]):
                continue
            n = normalize_edge(edge)
            idx = create_face(n)
            if idx is None:
                continue
            face = faces[idx]
            if is_original_point(edge[0]):
                face.add_neighbor(edge[0])
            if is_original_point(edge[1]):
                face.add_neighbor(edge[1])

    cells = [[] for i in range(len(orig_points))]
    for face_idx in range(len(faces)):
        face = faces[face_idx]
        for neighbor in face.neighbors:
            cells[neighbor].append(face_idx)

def print_results(detailed=False):
    for c in range(len(cells)):
        cell = cells[c]
        print("Cell %s at %s (%d faces)" % (get_name(orig_points[c]), orig_points[c], len(cell)))
        for face_idx in cell:
            face = faces[face_idx]
            neighbors = set(face.neighbors)
            neighbors.remove(c)
            neighbors = [get_name(orig_points[n]) for n in list(neighbors)]
            print("\tF%03d %s: " % (face_idx, neighbors), end='')
            if not detailed:
                for vertex_name in face.vertices:
                    print("%s " % vertex_name, end='')
            else:
                print()
                for vertex_name in face.vertices:
                    print("\t\t%s %s" % (vertex_name, svec(name_to_vec[vertex_name])))
            print()
        print()

def remove_cocentric():
    to_remove = []
    for i in range(len(tetrahedra)-1):
        for j in range(i+1, len(tetrahedra)):
            diff = centers[i] - centers[j]
            norm = linalg.norm(diff)
            if norm < 1e-5:
                to_remove.append(j)

    to_remove.sort()
    for i in reversed(to_remove):
        del tetrahedra[i]
        del centers[i]

def print_vectors():
    print()

    for vec,name in names.items():
        print("%3s\t%s" % (name, vec))

def main(folder, name):
    global all_points, tetrahedra, neighbors, centers, orig_points, edges, faces
    all_points, tetrahedra, neighbors = read_tetgen(folder, name)
    num_orig = int((len(all_points) - 4) / 7)
    orig_points = all_points[:num_orig]
    save_point_names()
    centers = calc_centers()
    edges = build_edges()
    voronoi()
    print_results()
#    print_vectors()

if __name__=='__main__':
    main('points', 'all')
    #v = np.array(name_to_vec["V47"])
    #u = np.array(name_to_vec["V70"])
    #diff = v-u
    #norm = linalg.norm(diff)
    #print(norm)