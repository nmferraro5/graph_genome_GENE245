#!/usr/bin/env python

import csv
out_rows = []

vcf_file = 'ALL.chr2.sorted.vcf'
with open(vcf_file, 'r') as iF:
  lines = csv.reader(iF, delimiter='\t')
  for line in lines:
    start = int(line[1])
    end = start + 1
    new_row = [line[0], line[1], end, line[2]]
    out_rows.append(new_row)

with open('ALL.chr2.sorted.position.vcf', 'w') as oF:
  writer = csv.writer(oF, delimiter='\t')
  for row in out_rows:
    writer.writerow(row)
