#!/usr/bin/python

import itertools as it

def quote(content):
  if content is None:
    return ""
  s = str(content)
  return "'{0}'".format(s.replace("'","''"))

class Edge:
  def __init__(self, dest, length):
    self.dest = dest
    self.length = length

class Node:
  def __init__(self, content, *edges):
    self.content = content
    self.edges = edges
  
  def content(self):
    return self.content 

  def children(self):
    return [edge.dest for edge in self.edges]

  def lengths(self):
    return [edge.length for edge in self.edges]

  def __str__(self):
    s = ['{subtree}{length}'.format(
      subtree=edge.dest,
      length=(":" + str(edge.length)) if edge.length is not None else '')
      for edge in self.edges]
    s = '(' + ','.join(s) + ')'
    return s + quote(self.content)

class Leaf(Node):
  def __init__(self, content):
    self.content = content
    self.edges = []

  def __str__(self):
    return quote(self.content)

# TODO(brunokim): Implement this without __str__ to avoid stack overflow
def newick_format(tree):
  return str(tree) + ';'

def leafs(tree):
  ls = []
  stack = [tree]
  
  while len(stack) > 0:
    t = stack.pop()
    for child in t.children():
      if len(child.children()) > 0:
        stack.append(child)
      else:
        ls.append(child)

  return ls

try:
  from igraph import Graph
  def to_graph(tree, g=None, serial=None, is_unrooted=True):
    is_root = False
    if g is None:
      is_root = True
      g = Graph()

    if serial is None:
      serial = it.count(-1, -1)

    c = tree.content if tree.content is not None else str(serial.next())

    if is_root and is_unrooted and len(tree.children()) == 2:
      left, right = [to_graph(child, g, serial) for child in tree.children()]
      length = sum(tree.lengths())
      g.add_edge(left, right, length=length)
    else:
      g.add_vertex(name = c)
      for child, length in zip(tree.children(), tree.lengths()):
        child_content = to_graph(child, g, serial)
        g.add_edge(c, child_content, length=length)

    if is_root:
      return g
    return c

  def distance_matrix(tree):
    g = to_graph(tree)
    ls = g.es["length"]
    vs = [leaf.content for leaf in leafs(tree)]

    m = []
    for v in vs:
      epaths = g.get_shortest_paths(v, vs, weights=ls, output="epath")
      distances = [
          sum(ls[edge_index] for edge_index in epath)
        for epath in epaths
      ]
      m.append(distances)

    return m, vs

except ImportError:
  def to_graph(tree):
    raise ImportError('igraph is not installed')

  def distance_matrix(tree):
    raise NotImplementedError(
        'Current implementation needs igraph, which is not installed')

def test_tree():
  """
              D
              |
              |3
              |    4
             (z)- - - - E
  A           |
  |1    5     |2      7
 (x)- - - - -(y)- - - - - - - C
  |
  |3
  |
  B"""
  e_Ax = Edge(Leaf('A'), 1)
  e_Bx = Edge(Leaf('B'), 3)
  x = Node('x', e_Ax, e_Bx)

  e_Dz = Edge(Leaf('D'), 3)
  e_Ez = Edge(Leaf('E'), 4)
  z = Node('z', e_Dz, e_Ez)

  e_xy = Edge(x, 5)
  e_zy = Edge(z, 2)
  e_Cy = Edge(Leaf('C'), 7)
  y = Node('y', e_xy, e_zy, e_Cy)
  
  m = [
      [ 0,  4, 13, 11, 12],
      [ 4,  0, 15, 13, 14],
      [13, 15,  0, 12, 13],
      [11, 13, 12,  0,  7],
      [12, 14, 13,  7,  0]
  ]

  return y, m, ['A', 'B', 'C', 'D', 'E']

if __name__ == '__main__':
  t, expected_m, expected_ids = test_tree()
  ls = leafs(t)
  ids = [l.content for l in ls]
  assert set(expected_ids) == set(ids)
  permutation = [ids.index(expected) for expected in expected_ids]

  m, ids = distance_matrix(t)
  n = len(ids)
  m = [[m[i][j] for j in permutation] for i in permutation]
  assert m == expected_m
