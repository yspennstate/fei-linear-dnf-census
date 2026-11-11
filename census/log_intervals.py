"""Certified dyadic enclosures of log2(integer), using positive atanh series.

Every operation is integer arithmetic. Series tail and every rounded product
are enclosed. Outputs floor/ceiling in units of 2^-52 bits.
"""
import json
import time
P=1<<144
S=1<<52
TERMS=52

def divceil(a,b):
    return (a+b-1)//b

def ln_ratio(num,den):
    # ln(num/den)=2 atanh((num-den)/(num+den)); 1<=num/den<=2.
    zlo=(num-den)*P//(num+den)
    zhi=divceil((num-den)*P,num+den)
    z2lo=zlo*zlo//P
    z2hi=divceil(zhi*zhi,P)
    plo,phi=zlo,zhi
    lo=hi=0
    for j in range(TERMS):
        lo+=plo//(2*j+1)
        hi+=divceil(phi,2*j+1)
        plo=plo*z2lo//P
        phi=divceil(phi*z2hi,P)
    # Now phi encloses P*z^(2TERMS+1). Geometric tail with all
    # denominators replaced by 2TERMS+1 is an upper bound.
    hi+=divceil(phi*P,(2*TERMS+1)*(P-z2hi))
    return 2*lo,2*hi

def main():
    start=time.time()
    l2lo,l2hi=ln_ratio(2,1)
    with open('log2_intervals.txt','w') as f:
        f.write(f'{S}\n0 0\n')
        for r in range(1,65537):
            k=r.bit_length()-1
            lo,hi=ln_ratio(r,1<<k)
            a=k*S+lo*S//l2hi
            b=k*S+divceil(hi*S,l2lo)
            assert 0<=b-a<=1
            f.write(f'{a} {b}\n')
    json.dump({'scale':S,'internal_bits':144,'terms':TERMS,
               'max_log_argument':65536,'seconds':time.time()-start},
              open('log_intervals.json','w'),indent=2)
if __name__=='__main__':
    main()
