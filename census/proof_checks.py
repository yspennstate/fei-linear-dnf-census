import collections
import itertools
import json
import math
import os
import random
import time
import numpy as np
from analyze import wht

assert os.getpriority(os.PRIO_PROCESS,0)>=10

def eta(a):
    a=np.asarray(a,dtype=float)
    out=np.zeros_like(a)
    z=a>0
    out[z]=-a[z]*np.log2(a[z])
    return out

def measure(n,ts,coords=None):
    x=np.arange(1<<n,dtype=np.int64)
    tt=np.zeros(1<<n,dtype=bool)
    for t in ts:tt|=(x&t)==t
    a=wht(1-2*tt.astype(np.int64))
    mu=(a/(1<<n))**2
    coords=range(n) if coords is None else coords
    er={}
    for i in coords:
        lo=x[(x&(1<<i))==0]
        p,q=mu[lo],mu[lo|(1<<i)]
        er[i]=(float(q.sum()),float((eta(p)+eta(q)-eta(p+q)).sum()))
    return float(tt.mean()),float(eta(mu).sum()),er

def mix(ts,i,k):
    incident=[t for t in ts if t&(1<<i)]
    blocks=[[j for j in range(max(ts).bit_length()) if (t^(1<<i))&(1<<j)] for t in incident]
    if any(len(b)<k for b in blocks):return None
    aa=[2.**-len(b) for b in blocks];b=sum(aa)
    if b>.25:return None
    loads={j:sum(2.**(1-t.bit_count()) for t in ts if t&(1<<j) and not t&(1<<i)) for block in blocks for j in block}
    termloads=[sum(loads[j] for j in block) for block in blocks]
    theta=1-(sum(a*l for a,l in zip(aa,termloads))+sum(a*c for a,c in itertools.combinations(aa,2)))/b
    if theta<=0:return None
    eps=min(1.,sum(a*sum(math.comb(len(block),r) for r in range(k)) for a,block in zip(aa,blocks)))
    al=0.
    for block in blocks:
        for J in itertools.combinations(block,k):
            den=1-.5*sum(loads[j] for j in J)
            if den<=0:return None
            al+=2.**(k*(k-1)//2-2*k)*math.prod(loads[j] for j in J)/den
    kap=al/(b*theta)
    err=2*eps+8/math.log(2)*math.sqrt(kap)+eps*math.log2(math.e*(3+8/(eps*b*theta)))+b*math.log2(math.e*(1+b**-2))
    return err

def main():
    start=time.time();short=rel=mc=0;worst_rel=worst_short=worst_mix=-1e30
    cases=[tuple(map(int,line.split())) for line in open('census6.families')]
    rng=random.Random(20261110)
    # Reservoir sample without keeping the full census in Python memory.
    sample=[]
    for j,line in enumerate(open('census8.families')):
        if j<200:sample.append(tuple(map(int,line.split())))
        else:
            k=rng.randrange(j+1)
            if k<200:sample[k]=tuple(map(int,line.split()))
    for n,ts in [(6,t) for t in cases]+[(8,t) for t in sample]:
        p,H,er=measure(n,ts)
        v=4*p*(1-p)
        for i,(d,e) in er.items():
            rel+=1;worst_rel=max(worst_rel,e/d-(2*abs(1-2*p)+math.log2(math.e*v/d)))
            for k in range(1,4):
                error=mix(ts,i,k)
                if error is not None:
                    mc+=1;worst_mix=max(worst_mix,abs(e/d-2)-error)
        for t in ts:
            w=t.bit_count();ii=max((i for i in er if t&(1<<i)),key=lambda i:er[i][0]);d,e=er[ii]
            short+=1
            assert d+1e-14>=2*(1-p)/(2**w-1)
            worst_short=max(worst_short,e/d-(2+math.log2(2*math.e*(2**w-1))))
    local=[]
    # Width-four term, plus a width-sixteen term meeting it in one bit.
    # The four coordinates satisfy LOCAL, including the shared coordinate.
    ts=[15,1|sum(1<<j for j in range(4,19))]
    p,H,er=measure(19,ts,range(4))
    for i,(d,e) in er.items():
        lam=sum(2.**(1-t.bit_count()) for t in ts if t&(1<<i))
        assert lam<=(1+1/(64*3**2))*2.**-3
        assert e<6*d
        local.append({'coordinate':i,'load':lam,'Inf':d,'E_Inf':e/d})
    lerrs=[]
    for L in range(3,1000):
        a=2.**-L;u=1/(64*L*L);theta=1-L*u*a
        eps=(math.sqrt(a)+math.sqrt(u))**2/theta
        kap=L*u/(4*theta*(1-u*a/2));b=(1+u)*a
        # Avoid inverse underflow by log-domain evaluation at large L.
        entropy=eps*(math.log2(math.e)+(math.log2(8+3*eps*a*theta)-math.log2(eps)-math.log2(a)-math.log2(theta)))
        empty=b*(math.log2(math.e)+math.log2(1+b*b)-2*math.log2(b))
        bound=2*eps+8/math.log(2)*math.sqrt(kap)+entropy+empty
        lerrs.append((bound,L))
    td=[]
    for width,q in [(3,3),(4,5)]:
        terms=[sum(1<<(x*q+(aa*x+bb)%q) for x in range(width)) for aa in range(q) for bb in range(q)]
        row=json.load(open(f'TD_{width}_{q}.json'));p,H,er=measure(width*q,terms,[0]);d,e=er[0]
        assert abs(H-row['H'])<1e-10
        assert abs(e/d-row['coordinate_orbits'][0]['E_Inf'])<1e-10
        td.append({'width':width,'q':q,'entropy_difference':H-row['H'],'erasure_difference':e/d-row['coordinate_orbits'][0]['E_Inf']})
    assert max(worst_rel,worst_short,worst_mix)<1e-10
    json.dump({'REL_checks':rel,'SHORT_checks':short,'MIX_checks':mc,
               'REL_max_signed_violation':worst_rel,'SHORT_max_signed_violation':worst_short,'MIX_max_signed_violation':worst_mix,
               'LOCAL_example':local,'LERR_max_on_3_to_999':max(lerrs),
               'TD_independent_numpy':td,'seconds':time.time()-start},open('proof_checks.json','w'),indent=2)
if __name__=='__main__':main()
