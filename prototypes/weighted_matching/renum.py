#!/usr/bin/env python3
# Renumber the named variables wsip writes (x<u>_<t>, y<k>_<t>_<s>) to x1, x2, ..., which is
# what RoundingSat accepts, and rewrite the OPB header counts to match.

import sys,re
ids={}
def r(m):
    n=m.group(0)
    if n not in ids: ids[n]=len(ids)+1
    return 'x%d'%ids[n]
lines=open(sys.argv[1]).read().split('\n')
out=[re.sub(r'\b[xy][0-9]+_[0-9_]+\b',r,l) for l in lines[1:]]
nc=sum(1 for l in out if l.strip() and not l.startswith('min'))
print('* #variable= %d #constraint= %d'%(len(ids),nc)); print('\n'.join(out))
