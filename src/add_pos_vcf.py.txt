#!/usr/bin/env python

import csv
import sys

out_rows = []

vcf_file = sys.argv[1]
out_file = sys.argv[2]

with open(vcf_file, 'r') as iF:
  lines = csv.reader(iF, delimiter='\t')
  for line in lines:
    try:
      start = int(line[1])
      end = start + 1
      new_row = [line[0], line[1], end, line[2]]
      out_rows.append(new_row)
    except:
      continue

with open(out_file, 'w') as oF:
  writer = csv.writer(oF, delimiter='\t')
  for row in out_rows:
    writer.writerow(row)
