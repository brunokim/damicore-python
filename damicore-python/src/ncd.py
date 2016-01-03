#!/usr/bin/python

import argparse
import os
import sys
import multiprocessing as mp
from shutil import copyfileobj as copy
from subprocess import Popen, PIPE, call
from progress_bar import ProgressBar

def gzip_compression(fname, slowness = 6, **kwargs):
  """Compression using gzip executable.

  @param fname Name of file to compress
  @param slowness Tradeoff parameter: 1 is fastest, 9 is best compression
  @return Size of compressed file in bytes
  """
  process = Popen(['gzip', '-c', '-%d' % slowness, fname], stdout=PIPE)
  compressed_size = len(process.communicate()[0])

  return compressed_size

def bzip2_compression(fname, slowness = 6, **kwargs):
  """Compression using bzip2 executable.
  
  @param fname Name of file to compress
  @param slowness Tradeoff parameter: 1 is fastest, 9 is best compression
  @return Size of compressed file in bytes
  """
  process = Popen(['bzip2', '-c', '-%d' % slowness, fname], stdout=PIPE)
  compressed_size = len(process.communicate()[0])

  return compressed_size

def ppmd_compression(fname, model_order = 6,
    ppmd_tmp_dir = os.path.join('ppmd_tmp'), **kwargs):
  """Compression using ppmd executable.

  @param fname Name of file to compress
  @param model_order Maximum model order to use, from 2 to 16
  @param ppmd_tmp_dir Temporary directory to use for ppmd output. It should be
      different from the input file directory, as the compressed file will have
      the same name as the input.
  @return Size of compressed file in bytes
  """
  tmp_fname = os.path.join(ppmd_tmp_dir, os.path.basename(fname))
  print fname, tmp_fname
  with open(os.devnull, 'w') as devnull:
    call(['ppmd', 'e', '-o%d' % model_order, '-f%s' % tmp_fname, fname],
        stdout=devnull)
  compressed_size = os.path.getsize(tmp_fname)
  os.remove(tmp_fname)

  return compressed_size

## Available compression functions
compression = {
    'gzip': gzip_compression,
    'bzip2': bzip2_compression,
    'ppmd': ppmd_compression,
}

#### Pairing functions ####

def concat(fname1, fname2, pair_dir = os.path.join('tmp'), **kwargs):
  """Concatenates given files into a new file.

  @param fname1, fname2 Files to concatenate (in this order)
  @param pair_dir Directory to output paired file
  @return Name of concatenated file
  """
  concat_name = os.path.basename(fname1) + '++' + os.path.basename(fname2)
  concat_fname = os.path.join(pair_dir, concat_name)

  with open(concat_fname, 'wb') as f,\
      open(fname1, 'rb') as i1,\
      open(fname2, 'rb') as i2:
    copy(i1, f)
    copy(i2, f)
  return concat_fname

def interleave(fname1, fname2, block_size = 1024,
    pair_dir = os.path.join('tmp'), **kwargs):
  """Interleaves blocks of the given files into a new file.

    Partitions each file into blocks with the given size, except for the last
  block, which can be smaller. Then, the blocks are interleaved in turn to
  create a new file. If the files have different sizes, the remaining blocks
  from the bigger one are concatenated.

      f1 = [--x1--][--x2--][--x3--][--x4--][--x5-]
      f2 = [--y1--][--y2--][-y3-]
      f = [--x1--][--y1--][--x2--][--y2--][--x3--][-y3-][--x4--][--x5-]

      This procedure produces NCD closer to zero for NCD(x,x), if x is bigger
  than the size limit for the compressor (32 KiB for gzip, 900 KiB for bzip2).

  @param fname1, fname2 Files to interleave (in this order)
  @param block_size Block size (in bytes)
  @param pair_dir Directory to output paired file
  @return Name of paired file
  """
  inter_name = (os.path.basename(fname1)
      + '--' + os.path.basename(fname2)
      + '--' + str(block_size))
  fname = os.path.join(pair_dir, inter_name)

  maxsize = max(os.path.getsize(fname1), os.path.getsize(fname2))
  with open(fname, 'wb') as f,\
      open(fname1, 'rb') as i1,\
      open(fname2, 'rb') as i2:
    for _ in xrange(0, maxsize, block_size):
      x1 = i1.read(block_size)
      x2 = i2.read(block_size)
      f.write(x1 + x2)
  return fname

## Available pairing functions
pairing = {
    'concat': concat,
    'interleave': interleave,
}

#### NCD functions ####

class NcdResult:
  def __init__(self, x, y, zx, zy, zxy, ncd):
    self.x = x
    self.y = y
    self.zx = zx
    self.zy = zy
    self.zxy = zxy
    self.ncd = ncd

def ncd(compression_fn, pairing_fn, fname1, fname2, compressed_sizes = None,
    **kwargs):
  """NCD calculation for a given pair of files.

      The normalized compression distance (NCD) between a pair of objects (x, y) 
  is defined as

                   Z(p(x,y)) - min{ Z(x), Z(y) }
      NCD_Z(x,y) = -----------------------------
                        max{ Z(x), Z(y) }

      where Z is the size of the compression of a given object and p is a
  pairing function that creates an object from two others. Theoretically, this
  distance is normalized between 0 and 1:
      NCD(x,x) == 0 and
      NCD(x,y) == 1 <=> Z(x) + Z(y) == Z(p(x,y))

  @param compression_fn Compression function with type
      fname, **kwargs -> compression_size
  @param pairing_fn Pairing function with type
      fname1, fname2, **kwargs -> paired_fname
      with the side-effect of creating a paired file
  @param fname1, fname2 Names of files to compare
  @param compressed_sizes Optional 2-tuple containing the compressed sizes of
      fname1 and fname2. If not provided, the sizes are calculated using
      compression_fn
  @param kwargs Additional arguments that are passed to compression_fn and
      pairing_fn
  @return NcdResult object containing additional information about the
      calculation
  """
  if compressed_sizes is None:
    compressed_sizes = (
        compression_fn(fname1, **kwargs),
        compression_fn(fname2, **kwargs))

  fname = pairing_fn(fname1, fname2, **kwargs)
  paired_compressed_size = compression_fn(fname, **kwargs)
  os.remove(fname)

  minimum, maximum = sorted(compressed_sizes)
  result = NcdResult(
      x = os.path.basename(fname1), y = os.path.basename(fname2),
      zx = compressed_sizes[0], zy = compressed_sizes[1],
      zxy = paired_compressed_size,
      ncd = float(paired_compressed_size - minimum)/maximum)
  return result

#### Parallel wrappers ####

def _parallel_compression_worker(args):
  """Wrapper for parallel calculation of compressed sizes."""
  compression_name, fname, queue, progress_bar, kwargs = (
      args.get('cname'), args.get('fname'), args.get('queue'),
      args.get('progress'), args.get('kwargs'))

  if compression_name is None:
    raise Exception('Compression not given')
  if fname is None:
    raise Exception('Filename not given')

  compression_fn = compression[compression_name]
  x = compression_fn(fname, **kwargs)

  if queue is not None:
    queue.put(x)

  if progress_bar is not None:
    progress_bar.increment()

  return x

def _parallel_ncd_worker(args):
  """Wrapper for parallel calculation of NCD pairs."""
  compression_name, pairing_name, fname1, fname2, queue, progress_bar,\
  compressed_sizes, kwargs = (args.get('cname'), args.get('pname'),
      args.get('f1'), args.get('f2'),
      args.get('queue'), args.get('progress'), args.get('zip'),
      args.get('kwargs'))

  if compression_name is None:
    raise Exception('Compression function name not given')
  if pairing_name is None:
    raise Exception('Pairing function name not given')
  if fname1 is None or fname2 is None:
    raise Exception('Filenames not given')

  compression_fn = compression[compression_name]
  pairing_fn = pairing[pairing_name]
 
  result = ncd(compression_fn, pairing_fn, fname1, fname2, compressed_sizes,
      **kwargs)

  if queue is not None:
    queue.put(result)

  if progress_bar is not None:
    progress_bar.increment()

  return result

#### Distance matrix calculations ####

def _serial_distance_matrix(fnames, compression_fn, pairing_fn, **kwargs):
  """Serial calculation for distance matrix."""
  sys.stderr.write('Compressing individual files...\n')
  progress_bar = ProgressBar(len(fnames))
  
  def update_progress(fname):
    x = compression_fn(fname)
    progress_bar.increment()
    return x
  compressed_sizes = map(update_progress, fnames)

  zip_size = dict(zip(fnames, compressed_sizes))

  sys.stderr.write('\nCompressing file pairs...\n')
  file_pairs = [(fname1, fname2)
      for fname1 in fnames
      for fname2 in fnames
      if fname1 < fname2]
  progress_bar = ProgressBar(len(file_pairs))

  def update_progress(pair):
    fname1, fname2 = pair
    ncd_result = ncd(compression_fn, pairing_fn, fname1, fname2,
        (zip_size[fname1], zip_size[fname2]), **kwargs)
    progress_bar.increment()
    return ncd_result
  ncd_results = map(update_progress, file_pairs)

  sys.stderr.write('\n')
  return ncd_results

def _parallel_distance_matrix(fnames, compression_name, pairing_name, **kwargs):
  """Parallel calculation of distance matrix."""
  num_cpus = mp.cpu_count()
  manager = mp.Manager()
  pool = mp.Pool(num_cpus)
  queue = manager.Queue(2*num_cpus)

  sys.stderr.write('Compressing individual files...\n')
  progress_bar = ProgressBar(len(fnames))

  compression_args = [{
    'cname': compression_name, 'fname': fname,
    'queue': queue, 'kwargs': kwargs}
    for fname in fnames]
 
  async_result = pool.map_async(_parallel_compression_worker, compression_args)

  for _ in xrange(len(fnames)):
    queue.get(timeout=5)
    progress_bar.increment()

  compressed_sizes = async_result.get()

  zip_size = dict(zip(fnames, compressed_sizes))

  sys.stderr.write('\nCompressing file pairs...\n')
  file_pairs = [(fname1, fname2)
      for fname1 in fnames
      for fname2 in fnames
      if fname1 < fname2]
  progress_bar = ProgressBar(len(file_pairs))

  ncd_args = [{
    'cname': compression_name,
    'pname': pairing_name,
    'f1': fname1, 'f2': fname2,
    'queue': queue, 'zip': (zip_size[fname1], zip_size[fname2]),
    'kwargs': kwargs} for fname1, fname2 in file_pairs]

  async_result = pool.map_async(_parallel_ncd_worker, ncd_args)
  pool.close()

  for _ in xrange(len(file_pairs)):
    queue.get(timeout=5)
    progress_bar.increment()

  ncd_results = async_result.get()
  sys.stderr.write('\n')
  return ncd_results

def distance_matrix(directory, compression_name, pairing_name,
    is_parallel=True, **kwargs):
  """Calculates matrix of distances between all files in a given directory.

  @param directory Directory with files to compare
  @param compression_name Name of compression function to use
  @param pairing_name Name of pairing function to use
  @param is_parallel Whether to perform computation in parallel
  @param kwargs Additional arguments for compression and pairing functions
  @return List of NcdResult objects from comparing all pairs of files
  """
  fnames = sorted(os.listdir(directory))
  fnames = [os.path.join(directory, fname)
      for fname in fnames
      if os.path.isfile(os.path.join(directory, fname))]

  if is_parallel:
    ncd_results = _parallel_distance_matrix(fnames, compression_name,
        pairing_name, **kwargs)
  else:
    ncd_results = _serial_distance_matrix(fnames, compression[compression_name],
        pairing[pairing_name], **kwargs)

  return ncd_results

#### Formatting functions ####

def csv_format(ncd_results, header=True):
  """Formats a list of NcdResult objects in CSV format.

  @param ncd_results List of NcdResult as returned by distance_matrix
  @param header Whether a header should be outputted to the CSV file
  @return String in CSV format
  """
  s = '' if not header else 'x,y,zx,zy,zxy,ncd\n'
  for result in ncd_results:
    s += '"{r.x}","{r.y}",{r.zx},{r.zy},{r.zxy},{r.ncd:.15f}\n'.format(
      r=result)
  return s

def to_matrix(ncd_results):
  """Converts a list of NcdResult objects to an n x n matrix.

  @param ncd_results List of NcdResult as returned by distance_matrix
  @return (m, ids) distance matrix with corresponding IDs
  """
  files = [r.x for r in ncd_results] + [r.y for r in ncd_results]
  ids = sorted(set(files))
  n = len(ids)

  m = [[0.0 for _ in xrange(n)] for _ in xrange(n)]

  for result in ncd_results:
    i, j = ids.index(result.x), ids.index(result.y)
    m[i][j] = m[j][i] = result.ncd

  return m, ids

def phylip_format(ncd_results, alternative_ids = None):
  """Formats a list of NcdResult objects in Phylip format.
   
      The Phylip format is used in phylogenetic software to store distance
  matrices between taxa. Each taxon name is limited to 10 chars, so the
  IDs used in NCD are truncated to satisfy this restriction. The format is as
  follows:

      <number-of-taxons>
      <taxon-name 1> <d(t1,t1)> <d(t1,t2)> ... <d(t1,tn)>
      <taxon-name 2> <d(t2,t1)> <d(t2,t2)> ... <d(t2,tn)>
      ...
      <taxon-name n> <d(tn,t1)> <d(tn,t2)> ... <d(tn,tn)>

  @param ncd_results List of NcdResult as returned by distance_matrix
  @param alternative_ids Optional IDs to use as taxon name. This might be
    necessary if the truncation of file names results in duplicates.
  @return String with matrix in Phylip format
  """
  m, ids = to_matrix(ncd_results)
  if alternative_ids is not None:
    ids = alternative_ids
  names = ['{name:<10.10}'.format(name=id_) for id_ in ids]
  # TODO(brunokim): Find conflicts and solve them

  s = '%d\n' % len(m)
  for name,row in zip(names, m):
    xs = ' '.join('%.15f' % dij for dij in row)
    s += name + ' ' + xs + '\n'

  return s

#### Command-line interface parser ####

def cli_parser():
  """Returns CLI parser for script.
  
  This may be useful for other scripts willing to call this one.
  """
  parser = argparse.ArgumentParser(
      description='Calculates NCD matrix between objects')
  parser.add_argument('directory',
      help='Directory containing files to compare')

  parser.add_argument('-c', '--compressor', choices=compression.keys(),
      default='ppmd', help='Compressor to use (default: ppmd)')
  parser.add_argument('-P', '--pairing', choices=pairing.keys(),
      default='concat', help='Pairing method to use (default: concat)')
  parser.add_argument('-o', '--output', help='output file (default: stdout)')
  parser.add_argument('-f', '--format', choices=['csv', 'phylip'],
      help='Choose matrix format (default: csv)')

  compressor_group = parser.add_argument_group('Compressor options', 
      'Options to control compressor behavior')
  compressor_group.add_argument('--slowness', '--gzip-slowness',
      '--bzip2-slowness', default=6, type=int,
      help='(gzip, bzip2) slowness of compression (1-9): ' + 
      '1 is faster, 9 is best compression')
  compressor_group.add_argument('--model-order', '--ppmd-model-order',
      default=6, type=int,
      help='(ppmd) model order (2-16): 2 is faster, 16 is best')
  compressor_group.add_argument('--memory', '--ppmd-memory',
      default=10, type=int, help='(ppmd) maximum memory, in MiB (1-256)')
  compressor_group.add_argument('--block-size', '--interleave-block-size',
      default=1024, type=int,
      help='(interleave) block size for interleaving, in bytes')

  misc_group = parser.add_argument_group('General options')
  
  is_serial = misc_group.add_mutually_exclusive_group()
  is_serial.add_argument('--serial', action='store_true',
      help='Compute compressions serially')
  is_serial.add_argument('--parallel', action='store_true',
      help='Compute compressions in parallel (default)')

  misc_group.add_argument('-v', '--verbose', action='count',
      help='Verbose output. Repeat to increase verbosity level (default: 1)',
      default = 1)
  misc_group.add_argument('--no-verbose', action='store_true',
      help='Turn verbosity off')
  misc_group.add_argument('-V', '--version', action='version', version='0.0.1')

  return parser

if __name__ == '__main__':
  parser = cli_parser()
  a = parser.parse_args()

  # TODO(brunokim): refactor code to use verbosity level, probably using a
  # logging library
  # verbose=0: no output to stderr
  # verbose=1: print progress bars
  # verbose=2: print files being compressed and a 'x/total' progress info
  verbose = 0 if a.no_verbose else a.verbose
  if verbose != 1:
    sys.stderr.write('Note: verbosity level not implemented yet\n')

  if not os.path.exists('tmp') or not os.path.isdir('tmp'):
    os.mkdir('tmp')
  if a.compressor == 'ppmd' and (
      not os.path.exists('ppmd_tmp') or not os.path.isdir('ppmd_tmp')):
    os.mkdir('ppmd_tmp')
 
  kwargs = {
      'pair_dir': 'tmp',
      'ppmd_tmp_dir': 'ppmd_tmp',
      'slowness': a.slowness,
      'model_order': a.model_order,
      'memory': a.memory,
      'block_size': a.block_size,
  }
  
  results = distance_matrix(a.directory, a.compressor, a.pairing,
      is_parallel = not a.serial, verbosity_level = verbose, **kwargs)

  if a.format == 'phylip':
    out = phylip_format(results)
  else:
    out = csv_format(results)
  
  if a.output is None:
    print out
  else:
    with open(a.output, 'wt') as f:
      f.write(out)

