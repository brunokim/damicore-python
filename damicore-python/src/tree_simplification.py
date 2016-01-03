#!/usr/bin/python

from tree import Node, Leaf, Edge
from math import log10, ceil
from random import Random

def num_digits(x):
  """Returns number of decimal digits in x.
  
  >>> num_digits(99)
  2
  >>> num_digits(100)
  2
  >>> num_digits(101)
  3
  """
  return int(ceil(log10(x)))

def artificial_ids(n):
  """Generates n artificial ids."""
  return ['_n{num:>0{max}}'.format(num=i, max=num_digits(n))
        for i in xrange(n)]

def matrix_argmin(m):
  """Returns indices of the minimum value in m."""
  indices = xrange(len(m))
  return min([(i,j) for i in indices for j in indices],
      key = lambda (i,j): m[i][j])

def calculate_q(m, sums):
  """Calculates matrix Q of the neighbor joining algorithm.
  
      The value Q_{ij} gives the decrease of total sum of edge lengths in the
  tree if nodes i and j are joined.
  """
  n = len(m)
  indices = xrange(n)
  q = [[0.0 for _ in indices] for _ in indices]

  for i, row in enumerate(m):
    for j, dij in enumerate(row):
      q[i][j] = (n - 2) * dij - sums[i] - sums[j]

  for i in indices:
    q[i][i] = float("inf")

  return q

def update_distance_matrix(m, sums, i, j):
  """Returns a new distance matrix with nodes i and j joined.

      The tuple (i,j) must be sorted. A new distance matrix is created with
  nodes i and j joined under a new node X. The node is placed at index i, and
  all its distances are updated.
  
  @return (new_m, di, dj) new distance matrix; distances i->X and j->X
  """
  n = len(m)
  dij, si, sj = m[i][j], sums[i], sums[j]

  di = (dij + (si - sj)/(n - 2))/2
  dj = (dij + (sj - si)/(n - 2))/2
  dk = [(m[i][k] + m[j][k] - dij)/2 for k in xrange(n)]

  new_m = [[dij for dij in row] for row in m]

  for k in xrange(n):
    new_m[i][k] = new_m[k][i] = dk[k]

  new_m.pop(j)
  for row in new_m:
    row.pop(j)

  return new_m, di, dj

def join_neighbors(tree, i, j, di, dj):
  r"""Returns a new tree with nodes i and j joined under a new node.

      [ (n_0)  ...  (n_i) ... (n_j)  ...  (n_N) ] --->

      [ (n_0)  ...     ()      ...  (n_{N-1}) ]
                   di /  \
                     /    \ dj
                   (n_i)   \
                         (n_j)
  """
  new_tree = [node for node in tree]
  new_tree[i] = Node(None,
      Edge(tree[i], di),
      Edge(tree[j], dj))
  new_tree.pop(j)
  return new_tree

def neighbor_joining(m, ids=None):
  """Neighbor Joining algorithm.
  
      Given a distance matrix, the algorithm seeks a tree that approximates the
  measured distances by greedily choosing to join the pair of elements that
  minimize the total sum of edge lengths.
  """
  n = len(m)
  if ids is None:
    ids = artificial_ids(n)

  # Turn m symmetric (and floating-point) if it's not already
  m = [ 
      [(m[i][j] + m[j][i])/2.0 for i in xrange(n)]
      for j in xrange(n)]

  tree = [Leaf(id_) for id_ in ids]

  for _ in xrange(n, 2, -1):
    # Find closest neighbors
    s = map(sum, m)
    q = calculate_q(m, s)

    # Join neighbors and update distance matrix
    i, j = sorted(matrix_argmin(q))
    m, di, dj = update_distance_matrix(m, s, i, j)
    tree = join_neighbors(tree, i, j, di, dj)

  d = m[0][1]
  return join_neighbors(tree, 0, 1, d/2, d/2)[0]

def _random_joining(ids):
  """Generates a tree by repeatedly joining elements randomly."""
  r = Random()
  tree = [Leaf(id_) for id_ in ids]
  for _ in xrange(len(ids) - 1):
    i, j = sorted(r.sample(xrange(len(tree)), 2))
    tree = join_neighbors(tree, i, j, r.random(), r.random())

  return tree[0]

if __name__ == '__main__':
  from tree import test_tree
  expected_tree, m, ids = test_tree()
  tree = neighbor_joining(m, ids)

  print tree, expected_tree

  # Random test
  from tree import distance_matrix
  expected_tree = _random_joining(map(str, xrange(10)))
  m, ids = distance_matrix(expected_tree)
  tree = neighbor_joining(m, ids)
  print tree, expected_tree
