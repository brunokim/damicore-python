#!/usr/bin/python

import sys

class ProgressBar:
  def __init__(self, end, length = 20):
    self.end = end
    self.length = length
    self.count = 0
    sys.stderr.write(str(self))
  
  def increment(self):
    percent = 100 * self.count / self.end
    next_percent = 100 * (self.count + 1) / self.end

    delta = float(self.length) / self.end
    size = int(round(delta * self.count))
    next_size = int(round(delta * (self.count + 1)))

    if next_size != size or percent != next_percent: 
      sys.stderr.write(str(self))

    self.count += 1

  def __str__(self):
      x = float(self.count) / self.end
      progress = int(round(x * self.length)) * '#'
      return '\r[{progress:<{max}}] {percent:3.0f} %'.format(
          progress=progress, max=self.length, percent = x*100)

if __name__ == '__main__':
  pass
