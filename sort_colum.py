import os
import fileinput
for line in fileinput.input('InputName', inplace=1):
    line = line.rstrip(os.linesep)
    print(' '.join(sorted(line.split())))
