import sys,collections
def rd(f):
    d=[]
    for l in open(f):
        l=l.rstrip('\n')
        if l.startswith('>'): d.append([l[1:],''])
        elif d: d[-1][1]+=l
    return d
a=rd(sys.argv[1]);b=rd(sys.argv[2])
ca=collections.Counter(map(tuple,a));cb=collections.Counter(map(tuple,b))
if ca==cb: print(f"EQUIV(n={len(a)})"); sys.exit(0)
na=set(x[0] for x in a); 
diff=[n for n in na if collections.Counter(s for m,s in a if m==n)!=collections.Counter(s for m,s in b if m==n)]
print(f"DIFF n={len(a)}/{len(b)} records_differing={len(diff)} e.g. {[d[:30] for d in diff[:3]]}"); sys.exit(1)
