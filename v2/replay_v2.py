"""Independent certificate replay; --full also regenerates optimization proposals.

Default uses the retained, source-bound proposals as untrusted inputs and
recomputes exact/interval certificates and independent counting dynamics.
"""
from pathlib import Path
import argparse,json,subprocess,sys,os,time
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--full',action='store_true');args=p.parse_args()
env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';started=time.perf_counter();checks=[]
jobs=[]
if args.full:jobs+=[['convex_compile.py','64'],['statistics.py']]
jobs += [['tree_exact.py'],['noise_geometry_audit.py'],['precision_cert.py'],
         ['primal_certificate.py','64'],['dual_certificate.py','64','128'],
         ['aggregate_convex.py'],['statistics.py','--compiled',str(ROOT/'primal_cert_results.json')],
         ['noise_smooth_review.py']]
for job in jobs:
    t=time.perf_counter();r=subprocess.run([sys.executable,str(ROOT/job[0]),*job[1:]],cwd=ROOT,env=env,capture_output=True,text=True)
    row={'command':['python',*job],'pass':r.returncode==0,'seconds':time.perf_counter()-t};checks.append(row)
    print(' '.join(job), 'PASS'if row['pass']else'FAIL',flush=True)
    if not row['pass']:
        print(r.stdout[-3000:],r.stderr[-4000:]);raise SystemExit(r.returncode)
cv=json.loads((ROOT/'convex_results.json').read_text());assert cv['status']=='PASS'
assert all(x['relative_gain_lower']>.14 for x in cv['summary'])
out={'schema':'cyclecraft-v2-replay-v1','status':'PASS','full_proposal_regeneration':args.full,'checks':checks,'elapsed_seconds':time.perf_counter()-started}
(ROOT/'replay_v2_receipt.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
