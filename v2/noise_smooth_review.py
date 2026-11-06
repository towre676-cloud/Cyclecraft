"""Independent retained-sample check of the smooth graph-noise floors.

The main statistics code obtains the floor with a traffic Laplacian inverse
and adaptive integration.  This audit instead uses the edge-disjoint cycle
resistance formula and a periodic trapezoid rule on retained samples.
"""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parent

def main():
    source=json.loads((ROOT/'statistics_results.json').read_text());reviews=[]
    assert source.get('all_checks_pass') is True,'Wait for the complete statistics run before reviewing.'
    assert set(source['models'])=={'cycle','star','chain'}
    for name,model in source['models'].items():
        for label,case in model['cases'].items():
            t=case['trajectory'];p=np.array(t['p']);x=np.exp(np.array(t['E']))*p
            q=np.array(t['q']);edges=case['edges'];times=np.array(t['phase'])*case['period']
            traffic=np.column_stack([q[:,e]*(x[:,i]+x[:,j])for e,(i,j)in enumerate(edges)])
            left=np.trapezoid(1/np.sum(1/traffic[:,:3],axis=1),times)
            right=np.trapezoid(1/np.sum(1/traffic[:,3:],axis=1),times)if name=='cycle'else 0.
            floor=np.diag([left,right]);V=np.array(case['covariance_asymptotic'])
            excess=np.linalg.eigvalsh((V+V.T)/2-floor)
            discrepancy=np.max(abs(floor-np.array(case['poisson_check']['graph_floor'])))
            reviews.append({'architecture':name,'case':label,'independent_cycle_floor':floor.tolist(),
                            'minimum_excess_eigenvalue':float(excess.min()),
                            'floor_quadrature_discrepancy':float(discrepancy),
                            'pass':bool(excess.min()>-1e-10 and discrepancy<1e-8)})
    out={'schema':'cyclecraft-v2-independent-smooth-floor-v1',
         'status':'PASS'if all(x['pass']for x in reviews)else'FAIL',
         'method':'Independent cycle-resistance formula and periodic trapezoid integration; no Laplacian inverse or main-code corrector.',
         'cases':len(reviews),'reviews':reviews}
    (ROOT/'noise_smooth_review.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({'status':out['status'],'cases':out['cases'],
                      'maximum_floor_discrepancy':max(x['floor_quadrature_discrepancy']for x in reviews)},indent=2))
    assert out['status']=='PASS'

if __name__=='__main__':main()
