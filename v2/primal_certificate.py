"""Exact repaired occupation protocols and continuous-time feasibility proofs.

The gate and occupation-potential knots are rational. A constant potential
repair closes the probability path exactly; a rational initial condition makes
every state's average occupation exactly 1/5. Probability is cubic on each
segment. Polynomial constraints use exact Bernstein coefficients and entropy
uses outward-rounded Arb interval integration. Binary64 fields are for replay
and plotting only; exact_source identifies the certified protocol.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json, sys, time
import numpy as np
import sympy as sp
from flint import arb, fmpq, ctx

ROOT = Path(__file__).resolve().parent
ctx.prec = 180
EDGES = {
    'star': [(0,1),(1,2),(2,0),(0,3),(0,4)],
    'chain': [(0,1),(1,2),(2,0),(0,3),(3,4)],
    'cycle': [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],
}
THETA = F(63,64)
CENTER_X = F(3,16)
EBOX = F(6,5)
RMAX = F(6)
PMIN = F(3,200)
SLEW_CAP = F(20)


def aq(x):
    if isinstance(x, F):
        return arb(fmpq(x.numerator, x.denominator))
    return arb(x)


def endpoint(x):
    m, e = x.man_exp()
    return F(int(m)) * F(2)**int(e)


def encode(x):
    lo, hi = x.lower(), x.upper()
    def dyad(v):
        m, e = v.man_exp()
        return {'mantissa':str(m), 'exponent2':int(e)}
    return {'ball':x.str(56), 'lower_dyadic':dyad(lo), 'upper_dyadic':dyad(hi),
            'lower_float':float(lo), 'upper_float':float(hi)}


def hull(lo, hi):
    return aq((lo+hi)/2) + arb(0, aq((hi-lo)/2))


def dyadic_float(x, bits=40):
    return F(round(F(float(x))*2**bits), 2**bits)


def bernstein_bounds(coefficients, lo=F(0), hi=F(1)):
    """Exact range enclosure of a power-basis polynomial on [lo,hi]."""
    d = len(coefficients)-1
    h = hi-lo
    local = [sum((coefficients[i]*comb(i,j)*lo**(i-j)*h**j
                  for i in range(j,d+1)), F(0)) for j in range(d+1)]
    bern = [sum((local[i]*F(comb(k,i),comb(d,i)) for i in range(k+1)), F(0))
            for k in range(d+1)]
    return min(bern), max(bern)


def poly_add(a, b, factor=F(1)):
    n=max(len(a),len(b));out=[F(0)]*n
    for i,v in enumerate(a):out[i]+=v
    for i,v in enumerate(b):out[i]+=factor*v
    return out


def poly_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def integral(coefficients):
    return sum((v/F(i+1) for i,v in enumerate(coefficients)),F(0))


def stage_data(edges, q0, q1, x0, x1, dt):
    flux=[];activity=[];pdot=[[F(0)]*3 for _ in range(5)]
    for e,(i,j)in enumerate(edges):
        qe=[q0[e],q1[e]-q0[e]]
        diff=[x0[i]-x0[j],x1[i]-x0[i]-x1[j]+x0[j]]
        summ=[x0[i]+x0[j],x1[i]-x0[i]+x1[j]-x0[j]]
        fe=poly_mul(qe,diff);flux.append(fe);activity.append(poly_mul(qe,summ))
        for k,v in enumerate(fe):pdot[i][k]-=v;pdot[j][k]+=v
    return pdot,flux,activity


def repair(raw):
    kind=raw['architecture'];edges=EDGES[kind];N=raw['N'];T=F(str(raw['period']));dt=T/N
    q=[[dyadic_float(v,10) for v in row] for row in raw['q_knots']]
    if kind=='cycle':
        for row in q:row[3:]=[F(683,1024),F(683,1024),F(682,1024)]
    x=[[THETA*dyadic_float(v)+(1-THETA)*CENTER_X for v in row] for row in raw['x_knots']]
    end=[F(0)]*5;Lbar=sp.zeros(5)
    for k in range(N):
        j=(k+1)%N;pdot,flux,act=stage_data(edges,q[k],q[j],x[k],x[j],dt)
        for i in range(5):end[i]+=dt*integral(pdot[i])
        for e,(i,jj)in enumerate(edges):
            w=dt*(q[k][e]+q[j][e])/2
            z=sp.Rational(w.numerator,w.denominator)
            Lbar[i,i]+=z;Lbar[jj,jj]+=z;Lbar[i,jj]-=z;Lbar[jj,i]-=z
    # Adding c to all potential knots changes the endpoint by -(integral L)c.
    rhs=sp.Matrix([sp.Rational(v.numerator,v.denominator) for v in end[1:]])
    cr=Lbar[1:,1:].inv()*rhs
    correction=[F(0)]+[F(int(v.p),int(v.q))for v in cr]
    x=[[v+correction[i]for i,v in enumerate(row)]for row in x]
    # p(t)=p0+delta(t). Pick p0 so mean p is exactly uniform.
    delta=[F(0)]*5;mean_delta=[F(0)]*5;stage_derivatives=[]
    for k in range(N):
        j=(k+1)%N;pdot,_,_=stage_data(edges,q[k],q[j],x[k],x[j],dt)
        stage_derivatives.append(pdot)
        for i in range(5):
            coeff=[delta[i]]+[dt*pdot[i][a]/F(a+1)for a in range(3)]
            mean_delta[i]+=dt*integral(coeff)
            delta[i]+=dt*integral(pdot[i])
    assert all(v==0 for v in delta)
    p0=[F(1,5)-v/T for v in mean_delta]
    assert sum(p0,F(0))==1
    current=p0[:];pcoeff=[];counts=[F(0)]*len(edges);activity=F(0);pmean=[F(0)]*5
    for k in range(N):
        j=(k+1)%N;pdot,flux,acts=stage_data(edges,q[k],q[j],x[k],x[j],dt)
        coeffs=[]
        for i in range(5):
            coeff=[current[i]]+[dt*pdot[i][a]/F(a+1)for a in range(3)]
            coeffs.append(coeff);pmean[i]+=dt*integral(coeff)
            current[i]+=dt*integral(pdot[i])
        pcoeff.append(coeffs)
        for e in range(len(edges)):
            counts[e]+=dt*integral(flux[e]);activity+=dt*integral(acts[e])
    assert current==p0
    assert all(v==T/5 for v in pmean)
    assert all(sum(pcoeff[k][i][a]for i in range(5))==(1 if a==0 else 0)
               for k in range(N)for a in range(4))
    source={'q_knots':[[str(v)for v in row]for row in q],
            'x_knots':[[str(v)for v in row]for row in x],
            'p0':[str(v)for v in p0],
            'p_coefficients':[[[str(v)for v in c]for c in row]for row in pcoeff]}
    protocol={'architecture':kind,'mode':'assisted_certified','N':N,'period':float(T),
              'edges':[list(e)for e in edges],
              'q_knots':[[float(v)for v in row]for row in q],
              'x_knots':[[float(v)for v in row]for row in x],
              'p_knots':[[float(c[0])for c in row]for row in pcoeff],
              'p_coefficients':[[[float(v)for v in c]for c in row]for row in pcoeff],
              'exact_source':source}
    return protocol,q,x,pcoeff,counts,activity,correction


def evaluate_protocol(sol, phase):
    """Numerical replay of the exact rational protocol for independent ODE/FCS."""
    if 'protocol'in sol:sol=sol['protocol']
    N=sol['N'];T=sol['period'];dt=T/N;z=(float(phase)%1)*N;k=min(int(z),N-1);s=z-k;kp=(k+1)%N
    q0=np.asarray(sol['q_knots'][k]);q1=np.asarray(sol['q_knots'][kp]);q=q0+(q1-q0)*s
    x0=np.asarray(sol['x_knots'][k]);x1=np.asarray(sol['x_knots'][kp]);x=x0+(x1-x0)*s
    pc=np.asarray(sol['p_coefficients'][k]);p=pc@np.array([1,s,s*s,s*s*s])
    pdot=(pc[:,1]+2*pc[:,2]*s+3*pc[:,3]*s*s)/dt
    xdot=(x1-x0)/dt;E=np.log(x/p);E_dot=xdot/x-pdot/p
    D=np.zeros((5,len(sol['edges'])))
    for e,(i,j)in enumerate(sol['edges']):D[i,e]=-1;D[j,e]=1
    C=-q[:,None]*D.T*(x/p)[None,:];K=D@C
    pi=p/x;pi/=pi.sum();current=C@p
    return {'p':p,'x':x,'q':q,'E':E,'E_dot':E_dot,'K':K,'C':C,'pi':pi,
            'current':current,'p_dot':pdot,'x_dot':xdot}


def certify(raw, subdivisions=256):
    st=time.perf_counter();sol,q,x,pcoeff,counts,activity,correction=repair(raw)
    N=sol['N'];T=F(str(sol['period']));dt=T/N;edges=EDGES[sol['architecture']]
    e_minus_upper=endpoint((-aq(EBOX)).exp().upper())
    e_plus_lower=endpoint(aq(EBOX).exp().lower())
    polynomial_checks={'probability_floor':True,'lower_well_box':True,'upper_well_box':True,'escape_rate_ceiling':True}
    margins={k:None for k in polynomial_checks};p_lower=None;x_lower=None
    slew_upper=arb(0);slew_parts_upper=arb(0);energy_upper=arb(0);escape_upper=arb(0)
    # Exact Bernstein certificate on each of 64 subsegments of every knot span.
    bern_subdivisions=64
    for k in range(N):
        kp=(k+1)%N;dpoly=[[F(0),F(0)]for _ in range(5)]
        for e,(i,j)in enumerate(edges):
            for node in [i,j]:dpoly[node][0]+=q[k][e];dpoly[node][1]+=q[kp][e]-q[k][e]
        for i in range(5):
            pc=pcoeff[k][i];xc=[x[k][i],x[kp][i]-x[k][i]]
            polys={
                'probability_floor':poly_add(pc,[-PMIN]),
                'lower_well_box':poly_add(xc,pc,-e_minus_upper),
                'upper_well_box':poly_add([v*e_plus_lower for v in pc],xc,-1),
                'escape_rate_ceiling':poly_add([RMAX*v for v in pc],poly_mul(dpoly[i],xc),-1)}
            for j in range(bern_subdivisions):
                lo=F(j,bern_subdivisions);hi=F(j+1,bern_subdivisions)
                for name,poly in polys.items():
                    lower,_=bernstein_bounds(poly,lo,hi)
                    polynomial_checks[name] &= lower>=0
                    if margins[name]is None or lower<margins[name]:margins[name]=lower
                pl,pu=bernstein_bounds(pc,lo,hi);xl,xu=bernstein_bounds(xc,lo,hi)
                p_lower=pl if p_lower is None else min(p_lower,pl)
                x_lower=xl if x_lower is None else min(x_lower,xl)
                pb=hull(pl,pu);xb=hull(xl,xu)
                ep=(xb/pb).log();energy_upper=max(energy_upper,abs(ep).upper())
                dp=[F(a)*pc[a]/dt for a in range(1,4)];dl,du=bernstein_bounds(dp,lo,hi)
                sx=aq(xc[1]/dt)/xb;sp=hull(dl,du)/pb
                slew=sx-sp
                slew_parts_upper=max(slew_parts_upper,(abs(sx).upper()+abs(sp).upper()).upper())
                slew_upper=max(slew_upper,abs(slew).upper())
                d0,d1=bernstein_bounds(dpoly[i],lo,hi)
                escape_upper=max(escape_upper,(hull(d0,d1)*xb/pb).upper())
    # Outward-rounded interval quadrature of entropy; x and q are affine.
    entropy=arb(0)
    for k in range(N):
        kp=(k+1)%N
        for j in range(subdivisions):
            lo=F(j,subdivisions);hi=F(j+1,subdivisions)
            xb=[];qb=[]
            for i in range(5):
                a=x[k][i]+(x[kp][i]-x[k][i])*lo
                b=x[k][i]+(x[kp][i]-x[k][i])*hi
                xb.append(hull(min(a,b),max(a,b)))
            for e in range(len(edges)):
                a=q[k][e]+(q[kp][e]-q[k][e])*lo
                b=q[k][e]+(q[kp][e]-q[k][e])*hi
                qb.append(hull(min(a,b),max(a,b)))
            f=sum((qb[e]*(xb[i]-xb[jj])*(xb[i]/xb[jj]).log()
                   for e,(i,jj)in enumerate(edges)),arb(0))
            entropy+=f*aq(dt/subdivisions)
    entropy_budget=F(str(raw['entropy_budget']));activity_budget=F(str(raw['activity_budget']))
    checks={**polynomial_checks,
            'entropy_upper_below_budget':entropy.upper()<=aq(entropy_budget),
            'exact_activity_below_budget':activity<=activity_budget,
            'finite_well_slew_below_cap':slew_upper<=aq(SLEW_CAP),
            'uniform_mixture_frontier_slew_below_cap':slew_parts_upper<=aq(SLEW_CAP),
            'exact_periodicity':True,'exact_probability_normalization':True,
            'exact_uniform_mean_inventory':True,
            'exact_active_edge_counts_agree':counts[0]==counts[1]==counts[2],
            'exact_zero_auxiliary_net_counts':all(v==0 for v in counts[3:])}
    numeric= evaluate_protocol(sol,.317)
    checks['replay_master_equation_identity']=float(abs(numeric['K']@numeric['p']-numeric['p_dot']).max())<1e-12
    out={'architecture':sol['architecture'],'input_N':N,'status':'PASS'if all(checks.values())else'FAIL',
         'source_proposal_objective':raw['objective'],
         'shrink_factor':str(THETA),'shrink_center_x':str(CENTER_X),
         'constant_repair':[str(v)for v in correction],
         'budgets':{'period':float(T),'entropy':float(entropy_budget),'activity':float(activity_budget),
                    'energy_absolute':float(EBOX),'escape_rate':float(RMAX),
                    'probability_floor':float(PMIN),'energy_slew':float(SLEW_CAP),
                    'mean_probability':[.2]*5},
         'objective_exact':str(counts[0]),'objective':float(counts[0]),
         'edge_counts_exact':[str(v)for v in counts],
         'activity_exact':str(activity),'activity':float(activity),
         'entropy_enclosure':encode(entropy),
         'minimum_probability_lower':str(p_lower),'minimum_x_lower':str(x_lower),
         'polynomial_constraint_minimum_margins':{k:str(v)for k,v in margins.items()},
         'absolute_energy_interval_upper':encode(energy_upper),
         'escape_rate_interval_upper':encode(escape_upper),
         'well_slew_interval_upper':encode(slew_upper),
         'uniform_mixture_slew_upper':encode(slew_parts_upper),
         'entropy_subintervals_per_segment':subdivisions,
         'bernstein_subintervals_per_segment':bern_subdivisions,
         'arb_precision_bits':ctx.prec,'protocol':sol,'checks':{k:bool(v)for k,v in checks.items()},
         'elapsed_seconds':time.perf_counter()-st}
    return out


def main():
    st=time.perf_counter();N=int(sys.argv[1])if len(sys.argv)>1 else 32
    source=ROOT/f'convex_raw_{N}.json';raw=json.loads(source.read_text());results=[]
    for r in raw:
        if r['mode']!='assisted':continue
        out=certify(r)
        if not out['checks']['entropy_upper_below_budget']:
            out=certify(r,1024)
        results.append(out)
        print(r['architecture'],out['status'],out['objective'],out['activity'],
              out['entropy_enclosure']['upper_float'],out['well_slew_interval_upper']['upper_float'],flush=True)
    out={'schema':'cyclecraft.v2.continuous-primal-certificate.v1','status':'PASS'if all(r['status']=='PASS'for r in results)else'FAIL',
         'source_file':source.name,'precision':'Exact rational source and polynomial identities; outward-rounded Arb transcendental enclosures.',
         'results':results,'elapsed_seconds':time.perf_counter()-st}
    (ROOT/'primal_cert_results.json').write_text(json.dumps(out,indent=2)+'\n')
    assert out['status']=='PASS',[(r['architecture'],r['checks'])for r in results]


if __name__=='__main__':main()
