import collections
import functools
import itertools
import json
import math
import os
import time
import numpy as np
import networkx as nx

assert os.getpriority(os.PRIO_PROCESS,0)>=10

def wht(f):
    a=np.asarray(f,dtype=np.int64).copy()
    h=1
    while h<len(a):
        b=a.reshape(-1,2,h)
        u=b[:,0,:].copy()
        v=b[:,1,:].copy()
        b[:,0,:]=u+v
        b[:,1,:]=u-v
        h*=2
    return a

def truth(n,terms):
    return sum(1<<x for x in range(1<<n) if any(x&t==t for t in terms))

def restrict(n,tt,i,b):
    return sum(((tt>>((x&((1<<i)-1))|((x>>i)<<(i+1))|(b<<i)))&1)<<x for x in range(1<<(n-1)))

def eta(x):
    return -x*math.log2(x) if x>0 else 0.

@functools.lru_cache(None)
def stats(n,tt):
    a=wht([1-2*((tt>>x)&1) for x in range(1<<n)])
    m=a*a
    den=1<<(2*n)
    mu=m/den
    H=sum(eta(p) for p in mu)
    inf=[]
    er=[]
    for i in range(n):
        e=d=0.
        for x in range(1<<n):
            if x&(1<<i):continue
            p,q=mu[x],mu[x|(1<<i)]
            e+=eta(p)+eta(q)-eta(p+q)
            d+=q
        inf.append(d);er.append(e)
    return H,tuple(inf),tuple(er),tuple(int(x) for x in a)

@functools.lru_cache(None)
def opt(n,tt):
    H,inf,er,_=stats(n,tt)
    if H==0 and not any(inf):return 0.
    if not n:return 0.
    return min(er[i]+.5*opt(n-1,restrict(n,tt,i,0))+.5*opt(n-1,restrict(n,tt,i,1))
               for i in range(n) if inf[i])

def absorb(ts):
    s=sorted(set(ts),key=lambda t:(t.bit_count(),t))
    out=[]
    for t in s:
        if not any(t&u==u for u in out):out.append(t)
    return tuple(out)

def term_restrict(ts,i,b):
    out=[]
    for t in ts:
        if t&(1<<i):
            if not b:continue
            t^=1<<i
        t=(t&((1<<i)-1))|((t>>(i+1))<<i)
        out.append(t)
    return absorb(out)

def term_tree(n,terms,early=False,high_degree=False):
    @functools.lru_cache(None)
    def run(n,ts,pending):
        if not ts or ts==(0,):return (0.,0.,0.)
        if not pending:
            t=min(ts,key=lambda x:(x.bit_count(),x))
            deg={j:sum(bool(s&(1<<j)) for s in ts) for j in range(n)}
            pending=tuple(sorted([j for j in range(n) if t&(1<<j)],
                                 key=lambda j:((-deg[j] if high_degree else deg[j]),j)))
        i=pending[0]
        children=[]
        tt=truth(n,ts)
        H,inf,er,_=stats(n,tt)
        for b in (0,1):
            sub=term_restrict(ts,i,b)
            pp=tuple(j-(j>i) for j in pending[1:])
            if early and b==0:pp=()
            children.append(run(n-1,sub,pp))
        H0=stats(n-1,restrict(n,tt,i,0))[0]
        H1=stats(n-1,restrict(n,tt,i,1))[0]
        J=H-.5*(H0+H1)-er[i]
        return tuple(v+.5*(children[0][k]+children[1][k]) for k,v in enumerate((er[i],J,inf[i])))
    return run(n,absorb(terms),())

def prime_expr(coeffs,den):
    out=collections.Counter()
    for x,c in coeffs.items():
        p=2
        while p*p<=x:
            while x%p==0:
                out[p]+=c;x//=p
            p+=1
        if x>1:out[x]+=c
    out={str(p):c for p,c in sorted(out.items()) if c}
    g=math.gcd(den,*out.values())
    return {'denominator':den//g,'log2_prime_coefficients':{p:c//g for p,c in out.items()}}

def details(n,terms,with_tree=True):
    tt=truth(n,terms)
    H,inf,er,a=stats(n,tt)
    I=sum(inf);m=[x*x for x in a];M=1<<(2*n)
    hcoeff=collections.Counter({M:sum(m)})
    for p in m:
        if p:hcoeff[p]-=p
    i=min((j for j in range(n) if inf[j]),key=lambda j:er[j]/inf[j])
    ecoeff=collections.Counter()
    for x in range(1<<n):
        if x&(1<<i):continue
        p,q=m[x],m[x|(1<<i)]
        for v,c in ((p+q,p+q),(p,-p),(q,-q)):
            if v:ecoeff[v]+=c
    graph=nx.Graph();graph.add_nodes_from(range(n))
    for t in terms:
        if t.bit_count()==2:graph.add_edge(*[i for i in range(n) if t&(1<<i)])
    out={'n':n,'terms':list(terms),'widths':sorted(t.bit_count() for t in terms),
         'degrees':sorted(sum(bool(t&(1<<j)) for t in terms) for j in range(n)),
         'P_true':tt.bit_count()/(1<<n),'mean_sign':1-2*tt.bit_count()/(1<<n),
         'H':H,'I':I,'H_I':H/I,'influences':inf,'E_Inf':[e/d if d else None for e,d in zip(er,inf)],
         'min_E_Inf':er[i]/inf[i],'H_I_exact':prime_expr(hcoeff,sum(x.bit_count()*m[x] for x in range(1<<n))),
         'min_E_Inf_exact':prime_expr(ecoeff,round(inf[i]*M)),
         'graph_components':[len(c) for c in nx.connected_components(graph)]}
    if with_tree:
        trees={}
        for early,high in itertools.product((False,True),repeat=2):
            e,j,d=term_tree(n,terms,early,high)
            assert abs(d-I)<1e-10 and abs(e+j-H)<1e-10
            trees[f'{"early" if early else "whole"}_{"high" if high else "low"}']={'E_I':e/I,'J_I':j/I}
        out['trees']=trees
        out['optimal_tree_E_I']=opt(n,tt)/I
    return out

def brute_labeled(n):
    # Completely labeled recursion over all allowed edges, with singleton
    # absorption checked directly. No isomorphism routine is used.
    candidates=list(range(1,1<<n));truths=set()
    def rec(pos,ts):
        truths.add(truth(n,ts))
        for p in range(pos,len(candidates)):
            t=candidates[p]
            if all((t&u).bit_count()<=1 and t&u!=t and t&u!=u for u in ts):
                rec(p+1,ts+(t,))
    rec(0,())
    full=(1<<(1<<n))-1
    live=[tt for tt in truths if tt not in (0,full) and all(restrict(n,tt,i,0)!=restrict(n,tt,i,1) for i in range(n))]
    return {'n':n,'labeled_including_dummy_and_false':len(truths),
            'essential_labeled':len(live),
            'max_H_I':max(stats(n,tt)[0]/sum(stats(n,tt)[1]) for tt in live),
            'max_min_E_Inf':max(min(e/d for e,d in zip(stats(n,tt)[2],stats(n,tt)[1])) for tt in live)}

def main():
    start=time.time();census=json.load(open('census8.json'));rows=[]
    for row in census['rows']:
        n=row['n'];seen=set()
        for objective in ('max_min_E_Inf','max_H_I'):
            for a in row[objective]['attainers']:
                terms=a['terms'];mask=0
                for t in terms:mask|=t
                live=[j for j in range(8) if mask&(1<<j)]
                terms=tuple(sorted(sum(1<<k for k,j in enumerate(live) if t&(1<<j)) for t in terms))
                if terms in seen:continue
                seen.add(terms);d=details(n,terms);d['objective']=objective;rows.append(d)
    # Direct Walsh coefficients via the dense character matrix independently
    # verify the transform on every census attainer.
    maxerr=0
    for d in rows:
        n=d['n'];tt=truth(n,d['terms'])
        f=[1-2*((tt>>x)&1) for x in range(1<<n)]
        slow=[sum(f[x]*(-1 if (x&s).bit_count()%2 else 1) for x in range(1<<n)) for s in range(1<<n)]
        assert tuple(slow)==stats(n,tt)[3]
        delta=max((e/dl-(2*abs(d['mean_sign'])+math.log2(math.e*(1-d['mean_sign']**2)/dl)))
                   for e,dl in zip(stats(n,tt)[2],stats(n,tt)[1]) if dl)
        maxerr=max(maxerr,delta)
    brute=[brute_labeled(n) for n in range(1,6)]
    json.dump({'attainers':rows,'independent_labeled_bruteforce':brute,
               'direct_character_transform_match':True,'REL_max_positive_violation':maxerr,
               'seconds':time.time()-start},open('analysis.json','w'),indent=2)
if __name__=='__main__':main()
