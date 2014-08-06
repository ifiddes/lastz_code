import os, sys, re
for line in sys.stdin:
    line = re.split(r'[ \t]', line.strip())
    if int(line[2]) - int(line[1]) > int(sys.argv[1]):
        print "\t".join(line)
