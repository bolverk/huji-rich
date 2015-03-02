__author__ = 'zmbq'
import numpy as np
from os import path

def _read_tetgen_file(filename, type=float):
    def parse_line(line, type):
        """ Parses a line into a list of numbers """
        if line[0] == '#':
            return None
        return np.array([type(s) for s in line.split()])

    """
    Reads a TetGen generated output file
    :param filename
    :return: A list of number lists - one per line, ignoring the first line
    """
    with open(filename, "r") as inf:
        lines = inf.readlines()
    lines = [l for l in lines if l[0]!='#' and l.strip()]  # Remove comments and blank lines
    numbers = [parse_line(l, type) for l in lines]
    count = int(numbers[0][0])
    if count!=len(numbers)-1:
        raise ValueError("File %s count is wrong" % filename)
    return numbers[1:]  # No need to include the count

def read_tetgen_points(filename):
    node_file = _read_tetgen_file(filename, float)
    points = []
    for line in node_file:
        points.append(line[1:4])
    return points

def read_tetgen(folder, lead_name):
    """
    Reads the TetGen Delaunay output
    :param folder: Folder name
    :param lead_name: Leader file name - lead.1.ele, lead.1.neigh, lead.1.node are all read
    :return: points, tetrahedra, neighbors
        points: array of points
        tetrahedra: array of tetrahedra (indices into points)
        neighbors: array with four neighbors for each tetrahedron (-1 if no neighbor)
    """
    filename = path.join(folder, lead_name) + ".1.";
    node_file = _read_tetgen_file(filename + "node", float)
    points = read_tetgen_points(node_file)

    tetrahedra_file = _read_tetgen_file(filename + "ele", int)
    tetrahedra = []
    for line in tetrahedra_file:
        tetrahedron = line[1:]   # First number is the index
        tetrahedron = [t-1 for t in tetrahedron]  # All indices as 1 based, adjust to 0
        tetrahedra.append(tetrahedron)

    neighbor_file = _read_tetgen_file(filename + "neigh", int)
    neighbors = []
    for line in neighbor_file:
        neighbor = line[1:]
        neighbor = [n-1 for n in neighbor]
        neighbors.append(neighbor)

    return points, tetrahedra, neighbors
