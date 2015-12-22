#!python

"""
  Merge several magres files into one.
"""
import sys
from magres.format import MagresFile

if __name__ == "__main__":
  magres_file_paths = sys.argv[1:]

  magres_files = [MagresFile(magres_file_path) for magres_file_path in magres_file_paths]

  merged_magres_file = MagresFile.merge(magres_files)

  print(str(merged_magres_file))
