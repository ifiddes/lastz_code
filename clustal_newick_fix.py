"""regex for fixing clustalomega guide trees"""

import sys, os, re
t = re.compile(r".+:\d+-.+:")
for line in sys.stdin:
    m = re.search(t, line)
    if m is not None:
        start, end = m.span()
        fpos = line.index(":")#find first colon; hacky
        line = list(line)
        line[fpos] = ";"
        line.insert(end-1, "'")
        line.insert(start, "'")

        sys.stdout.write("".join(line))
    else:
        sys.stdout.write(line)
    