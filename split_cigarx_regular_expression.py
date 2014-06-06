#!/usr/bin/env python2.7

import re
r = re.compile(r"[1-9]*[=XID]")

def parse_cigar(cigar):
    return re.findall(r,cigar)