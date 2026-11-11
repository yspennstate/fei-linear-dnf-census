import hashlib
import json
import math
import os
from pathlib import Path
import mpmath as mp

mp.mp.dps=90

def evaluate(expr):
    return sum(mp.mpf(c)*mp.log(int(p),2) for p,c in expr['log2_prime_coefficients'].items())/expr['denominator']

data=json.load(open('analysis.json'))
exact=[]
for d in data['attainers']:
    h=evaluate(d['H_I_exact']);e=evaluate(d['min_E_Inf_exact'])
    assert abs(h-d['H_I'])<mp.mpf('1e-12')
    assert abs(e-d['min_E_Inf'])<mp.mpf('1e-12')
    exact.append({'n':d['n'],'terms':d['terms'],'H_I':mp.nstr(h,85),'min_E_Inf':mp.nstr(e,85)})
json.dump(exact,open('exact_values.json','w'),indent=2)

checks=[]
for n in range(1,9):
    c=json.load(open(f'census{n}.certified.json'))
    l=json.load(open(f'census{n}.labeled_check.json'))
    assert c['count']==l['classes'] and l['match']
    row={'n':n,'count':c['count'],'count_match':True}
    for field in ['H_I','min_E_Inf']:
        cand=c[field]['candidates']
        if n>2:assert min(x['lower'] for x in cand)>c[field]['other_upper']
        vals=[x for x in exact if x['n']==n]
        value=max(mp.mpf(x[field]) for x in vals)
        for x in cand:
            # Exact attainer values are independently recovered from their
            # integer prime-log expressions, rather than the log table.
            assert mp.mpf(str(x['lower']))<=value<=mp.mpf(str(x['upper']))
        row[field]=mp.nstr(value,85)
    checks.append(row)

for x in json.load(open('graph_names.json')).values():assert x['isomorphic']
pc=json.load(open('proof_checks.json'))
assert max(pc[k] for k in ['REL_max_signed_violation','SHORT_max_signed_violation','MIX_max_signed_violation'])<0
codefiles=['census.cpp','certify.cpp','labeled_check.cpp','log_intervals.py','analyze.py','structured.cpp','proof_checks.py','PROOFS.md','LITERATURE.md']
hashes={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in codefiles}
total=sum(p.stat().st_size for p in Path('.').rglob('*') if p.is_file())
assert total<200_000_000
json.dump({'checks':checks,'sha256':hashes,'measured_folder_bytes':total,
           'priority':os.getpriority(os.PRIO_PROCESS,0),'C_width_2':2+math.log2(6*math.e)+2/math.log(2),
           'essential_classes_through_8':sum(x['count'] for x in checks),
           'classes_with_dummies_and_constants_on_8':2+sum(x['count'] for x in checks)},
          open('final_audit.json','w'),indent=2)
print('All final audits passed; exact_values.json and final_audit.json saved.')
