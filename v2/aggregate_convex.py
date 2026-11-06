"""Combine exact feasible witnesses with continuous dual upper certificates."""
from pathlib import Path
from fractions import Fraction as F
import json,math
from flint import arb,ctx
ctx.prec=150
ROOT=Path(__file__).resolve().parent
def aq(f):f=F(f);return arb(f.numerator)/arb(f.denominator)
def up(a):return math.nextafter(float(a.upper()),math.inf)
def down(a):return math.nextafter(float(a.lower()),-math.inf)
def dyadic(d):return F(d['mantissa'])*F(2)**int(d['exponent2'])

def main():
    p=json.loads((ROOT/'primal_cert_results.json').read_text());d=json.loads((ROOT/'dual_bounds_64_128.json').read_text());raw=json.loads((ROOT/'convex_raw_64.json').read_text())
    assert p['status']=='PASS'
    ds={x['architecture']:x for x in d};bs={x['architecture']:x for x in raw if x['mode']=='baseline'}
    summary=[];frontier=[]
    for x in p['results']:
        kind=x['architecture'];c=ds[kind];Q=F(x['objective_exact']);U=F(c['continuous_upper']);beta=F(c['entropy_multiplier'])
        entropyU=dyadic(x['entropy_enclosure']['upper_dyadic']);A=F(x['activity_exact']);qL=down(aq(Q));uU=up(aq(U));gap=aq(Q)-aq(U)
        assert gap>0
        summary.append({'architecture':kind,'period':8.,'entropy_budget':.15,'activity_budget':16.,'energy_bound':1.2,'escape_bound':6.,'well_slew_bound':20.,
            'mean_inventory':[.2]*5,'assisted_current_exact':str(Q),'assisted_current_lower':qL,'assisted_current':float(Q),
            'baseline_numerical_proposal':bs[kind]['objective'],'baseline_continuous_upper':uU,'certified_gap_lower':down(gap),
            'relative_gain_lower':down(aq(Q)/aq(U)-1),'entropy_upper':up(aq(entropyU)),'activity':float(A),
            'certificate_subintervals':c['subintervals'],'certificate_status':'PASS','primal_checks':x['checks'],
            'baseline_class':'All continuous positive occupation/potential controls with fixed gates, auxiliary E3=E4=0, bounds and mean inventory; no well slew restriction.'})
        for theta in [F(3,5),F(7,10),F(4,5),F(9,10),F(1)]:
            S=theta*F(3,20)
            # The original upper uses an entropy budget rounded outward from binary64 .15.
            Sbase=F(math.nextafter(.15,math.inf))
            uu=aq(U)+aq(beta)*(aq(S)-aq(Sbase))
            qq=aq(theta*Q);gg=qq-uu
            activity=theta*A+(1-theta)*15
            frontier.append({'architecture':kind,'theta_exact':str(theta),'theta':float(theta),'entropy_budget_exact':str(S),'entropy_budget':float(S),
                'activity_budget':16.,'assisted_current':float(theta*Q),'assisted_current_lower':down(qq),
                'baseline_continuous_upper':up(uu),'certified_gap':down(gg),'certified_gap_lower':down(gg),
                'relative_gain_lower':down(qq/uu-1),'assisted_entropy_upper':up(aq(theta*entropyU)),
                'assisted_activity':float(activity),'certificate_status':'strict domination'if gg>0 else'valid bounds; no domination claimed',
                'feasibility':'Exact convex mixture with uniform p=1/5,x=3/16; common well/rate boxes and activity by convexity; well slew by uniform mixture certificate.'})
    out={'schema':'cyclecraft.v2.bounded-continuous-synthesis.v1','status':'PASS',
        'scope':'Continuous-time fixed-gate global comparator bounds, exact rational assisted witnesses, rigorous Arb/LP residual certificates.',
        'budgets':{'period':'8','entropy':'3/20','activity':'16','well_absolute':'6/5','escape_rate':'6','probability_floor':'3/200','well_slew':'20','mean_occupation_each':'1/5'},
        'control_model':'Baseline all three active wells free; auxiliary wells fixed zero. Assistance releases both auxiliary wells. No fixed root potential or tied gate-well constraint is imposed in this v2 benchmark.',
        'installation':'Three active edges, auxiliary conductance sum exactly 2 in every architecture; all schedules continuous piecewise linear; only active edge01 is gated.',
        'summary':summary,'frontier_rows':frontier,'primal_source':'primal_cert_results.json','dual_source':'dual_bounds_64_128.json',
        'checks':{'all_primal_certificates_pass':all(x['status']=='PASS'for x in p['results']),
                  'all_strict_continuous_comparisons':all(x['certified_gap_lower']>0 for x in summary),
                  'all_mixture_slew_bounds':all(x['uniform_mixture_slew_upper']['upper_float']<20 for x in p['results']),
                  'all_frontier_entropy_bounds':all(x['assisted_entropy_upper']<x['entropy_budget']for x in frontier),
                  'all_frontier_activity_bounds':all(x['assisted_activity']<16 for x in frontier)}}
    assert all(out['checks'].values());(ROOT/'convex_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out['summary'],indent=2))
if __name__=='__main__':main()
