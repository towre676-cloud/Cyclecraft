"""Cyclecraft II: rigorous finite-horizon certificates for dyadic hold protocols.

All source q, exp(E), durations and tilt z=exp(chi) are exact rationals.
Arb matrix exponentials enclose the exact propagators. Periodic initialization
is certified by a Doeblin contraction/residual estimate, independently checked
against 100-decimal-digit mpmath exponentials.
"""
from pathlib import Path
from fractions import Fraction as F
import json,time,sys
import numpy as np
import mpmath as mp
from flint import arb,arb_mat,fmpq,ctx

ROOT=Path(__file__).resolve().parent
ctx.prec=320
N=5
CORE=[(0,1),(1,2),(2,0)]
GRAPHS={'cycle':CORE+[(0,3),(3,4),(4,0)],
        'star':CORE+[(0,3),(0,4)]}
TAU=F(2)
W=[[F(1,4),F(4),F(4),F(1,2),F(2)],
   [F(4),F(1,4),F(4),F(1),F(1)],
   [F(4),F(4),F(1,4),F(2),F(1,2)]]
QC=[[F(1,8),F(1,8),F(2)],
    [F(2),F(1,8),F(1,8)],
    [F(1,8),F(2),F(1,8)]]

def A(x):
    if isinstance(x,F):return arb(fmpq(x.numerator,x.denominator))
    return arb(x)

def dyad(x):
    """Serialize an exact dyadic Arb endpoint without decimal rounding."""
    m,e=x.man_exp();return {'mantissa':str(m),'exponent2':int(e)}

def enc(x):
    lo=x.lower();hi=x.upper()
    return {'ball':x.str(78),'lower_dyadic':dyad(lo),'upper_dyadic':dyad(hi),
            'lower_float':float(lo),'upper_float':float(hi),
            'radius_float':float(x.rad())}

def maximum(xs):return max(xs,key=lambda x:x.upper())
def minimum(xs):return min(xs,key=lambda x:x.lower())

def generator(graph,stage,z=F(1),edge=0,lift=None):
    edges=GRAPHS[graph];auxq=F(1,2)if graph=='cycle'else F(3,4)
    q=QC[stage]+[auxq]*(len(edges)-3)
    size=N if lift is None else N*(2*lift+1)
    K=[[F(0)for j in range(size)]for i in range(size)]
    levels=[0]if lift is None else range(-lift,lift+1)
    for k in levels:
        for e,(i,j)in enumerate(edges):
            for source,target,sign in [(i,j,1),(j,i,-1)]:
                rate=q[e]*W[stage][source]
                col=source if lift is None else (k+lift)*N+source
                K[col][col]-=rate
                inc=sign if e==edge else 0
                kk=k+inc
                if lift is not None and not -lift<=kk<=lift:continue
                row=target if lift is None else (kk+lift)*N+target
                K[row][col]+=rate*(z**inc)
    return K

def mat(rows):return arb_mat([[A(v)for v in row]for row in rows])
def monodromy(graph,z=F(1),edge=0,lift=None):
    size=N if lift is None else N*(2*lift+1)
    M=arb_mat([[int(i==j)for j in range(size)]for i in range(size)])
    for s in range(3):M=(mat(generator(graph,s,z,edge,lift))*A(TAU)).exp()*M
    return M

def periodic_certificate(M):
    Id=arb_mat([[int(i==j)for j in range(N)]for i in range(N)])
    B=M-Id
    for j in range(N):B[N-1,j]=1
    rhs=arb_mat(N,1,[0]*(N-1)+[1]);proposal=B.solve(rhs)
    # Exact dyadic normalized probability vector; last coordinate is exact.
    bits=ctx.prec-32
    pp=[]
    for i in range(N-1):
        m,e=proposal[i,0].mid().man_exp();f=F(int(m))*F(2)**int(e)
        pp.append(F(round(f*(2**bits)),2**bits))
    pp.append(F(1)-sum(pp,F(0)))
    p=[A(v)for v in pp]
    assert all(v.is_exact()and v>0 for v in p)
    pv=arb_mat(N,1,p);res=(M-Id)*pv
    rnorm=sum((abs(res[i,0])for i in range(N)),arb(0)).upper()
    m=minimum([M[i,j].lower()for i in range(N)for j in range(N)])
    assert m>0
    contraction=1-N*m
    delta=(rnorm/(N*m)).upper()
    checkmass=maximum([abs(sum((M[i,j]for i in range(N)),arb(0))-1)for j in range(N)])
    assert checkmass.contains(0)
    return pv,delta,{'p_proposal':[enc(v)for v in p],
        'min_entry_lower':enc(m),'doeblin_contraction_upper':enc(contraction),
        'stationary_residual_l1_upper':enc(rnorm),'initial_l1_error_upper':enc(delta),
        'column_sum_residual_enclosure':enc(checkmass)}

def mgf(Mz,n,p,delta):
    Mn=Mz**n
    row=[sum((Mn[i,j]for i in range(N)),arb(0))for j in range(N)]
    center=sum((row[j]*p[j,0]for j in range(N)),arb(0))
    err=(maximum([abs(v).upper()for v in row])*delta).upper()
    return center+arb(0,err)

def survival(Msafe,n,p,delta,b):
    size=N*(2*b+1);pv=arb_mat(size,1,[p[i-b*N,0]if b*N<=i<(b+1)*N else 0 for i in range(size)])
    out=Msafe**n*pv;v=sum((out[i,0]for i in range(size)),arb(0))
    # Substochastic propagation has l1 norm at most one.
    return v+arb(0,delta)

def mp_matrix(rows):return mp.matrix([[mp.mpf(v.numerator)/v.denominator for v in row]for row in rows])
def mp_monodromy(graph,z=F(1),edge=0,lift=None):
    size=N if lift is None else N*(2*lift+1);M=mp.eye(size)
    for s in range(3):M=mp.expm(mp_matrix(generator(graph,s,z,edge,lift))*mp.mpf(2))*M
    return M

def shortcut_monodromy(eps,killed=False):
    """Star with rare extra 3--4 channel; kill its jumps if requested."""
    M=arb_mat([[int(i==j)for j in range(N)]for i in range(N)])
    for s in range(3):
        K=generator('star',s)
        for source,target in [(3,4),(4,3)]:
            rate=eps*W[s][source];K[source][source]-=rate
            if not killed:K[target][source]+=rate
        M=(mat(K)*A(TAU)).exp()*M
    return M

def moment_monodromy(graph,edge):
    """Analytic first/second counting derivatives, no finite differences."""
    size=3*N;M=arb_mat([[int(i==j)for j in range(size)]for i in range(size)])
    u,v=GRAPHS[graph][edge]
    for stage in range(3):
        K=generator(graph,stage);B=[[F(0)for j in range(size)]for i in range(size)]
        for block in range(3):
            for i in range(N):
                for j in range(N):B[block*N+i][block*N+j]=K[i][j]
        for source,target,sgn in [(u,v,1),(v,u,-1)]:
            rate=K[target][source]
            B[N+target][source]=sgn*rate
            B[2*N+target][source]=rate
            B[2*N+target][N+source]=2*sgn*rate
        M=(mat(B)*A(TAU)).exp()*M
    return M

def finite_moments(M,n,p,delta):
    ini=arb_mat(3*N,1,[p[i,0]+arb(0,delta)if i<N else 0 for i in range(3*N)])
    result=M**n*ini
    mean=sum((result[N+i,0]for i in range(N)),arb(0))
    second=sum((result[2*N+i,0]for i in range(N)),arb(0))
    variance=second-mean**2
    return mean,variance

def main():
    st=time.perf_counter();mp.mp.dps=110;checks=[];out={}
    def check(name,passed,detail=None):checks.append({'name':name,'pass':bool(passed),'detail':detail})
    for graph in GRAPHS:
        M=monodromy(graph);p,delta,pc=periodic_certificate(M)
        MM=mp_monodromy(graph);BB=MM-mp.eye(N)
        for j in range(N):BB[N-1,j]=1
        pp=mp.lu_solve(BB,mp.matrix([0]*(N-1)+[1]))
        mp_error=max(abs(mp.mpf(p[i,0].str(105,radius=False))-pp[i])for i in range(N))
        check(graph+' periodic reference enclosed',mp_error<mp.mpf(delta.upper().str(105,radius=False)),mp.nstr(mp_error,35))
        check(graph+' Doeblin contraction',float(pc['doeblin_contraction_upper']['upper_float'])<1)
        tilts={};tail_candidates=[];aux_candidates=[]
        for z in [F(1,2),F(2,3),F(3,4),F(7,8)]:
            Mz=monodromy(graph,z,0);g=mgf(Mz,10,p,delta)
            tilts[str(z)]=enc(g);tail_candidates.append((z,g.upper()))
        bz,bg=min(tail_candidates,key=lambda v:float(v[1]))
        # At integer-valued N_L, P[N_L<=0] <= E[z^N_L], z<1.
        tail={'event':'target signed count <= 0 after 10 periods','tilt_z':str(bz),
              'certified_probability_upper':enc(minimum([arb(1),bg])),
              'independent_copies_10_aggregate_wrong_direction_upper':enc(minimum([arb(1),(bg**10).upper()]))}
        for z in [F(9,8),F(5,4),F(3,2),F(2)]:
            gp=mgf(monodromy(graph,z,3),10,p,delta)
            gm=mgf(monodromy(graph,1/z,3),10,p,delta)
            # Union bound for integer signed auxiliary count magnitude>=12.
            bnd=(gp+gm)/A(z**12)
            aux_candidates.append((z,bnd.upper()))
        az,ag=min(aux_candidates,key=lambda v:float(v[1]))
        aux_tail={'event':'absolute auxiliary edge count >=12 after 10 periods','tilt_z':str(az),
                  'certified_probability_upper':enc(minimum([arb(1),ag]))}
        if graph=='star':aux_tail['exact_structural_probability']=0
        safe={}
        for b in [1,2]:
            Ms=monodromy(graph,edge=3,lift=b);S=survival(Ms,10,p,delta,b)
            safe[str(b)]={'survival_probability':enc(S),'failure_probability':enc(1-S),
                          'exact_survival_by_bridge_theorem':graph=='star'}
            if graph=='star':check('star bridge survival b'+str(b),S.contains(1))
            else:check('cycle positive crossing risk b'+str(b),(1-S)>0)
        # Reference the genuinely tilted propagator with a nontrivial rational tilt.
        zref=F(3,4);g=mgf(monodromy(graph,zref,0),10,p,delta)
        mref=(mp_monodromy(graph,zref,0)**10)*pp;gref=sum(mref)
        refball=arb(mp.nstr(gref,109))
        check(graph+' tilted reference enclosed',g.contains(refball),mp.nstr(gref,75))
        moments={}
        for name,edge in [('target',0),('auxiliary',3)]:
            Mm=moment_monodromy(graph,edge);rows=[]
            for periods in [1,10,100]:
                mean,var=finite_moments(Mm,periods,p,delta)
                check(graph+' '+name+' variance n'+str(periods),var>0)
                if name=='auxiliary':check(graph+' mean fence n'+str(periods),mean.contains(0))
                if name=='auxiliary'and graph=='star':check('bridge variance bounded n'+str(periods),var<1)
                rows.append({'periods':periods,'mean':enc(mean),'variance':enc(var)})
            moments[name]=rows
        out[graph]={'periodic':pc,'target_lower_tail':tail,'target_tilted_mgfs':tilts,
                    'auxiliary_terminal_tail':aux_tail,'auxiliary_survival':safe,
                    'finite_horizon_moments':moments,
                    'independent_110_digit_mgf_reference':mp.nstr(gref,100),
                    'tilted_reference_ball':enc(g)}
    shortcuts=[]
    for eps in [F(1,1024),F(1,256),F(1,64)]:
        M=shortcut_monodromy(eps);p,delta,pc=periodic_certificate(M)
        Ms=shortcut_monodromy(eps,True)
        count=Ms**10*p;S=sum((count[i,0]for i in range(N)),arb(0))+arb(0,delta)
        risk=1-S
        integrated_cap=10*TAU*eps*sum(max(row[3],row[4])for row in W)
        cap=1-(-A(integrated_cap)).exp()
        check('rare shortcut hazard '+str(eps),risk>0 and cap-risk>0)
        shortcuts.append({'extra_edge':[3,4],'extra_conductance':str(eps),
            'initialization':'unique periodic law of the full leaky-star protocol',
            'any_shortcut_jump_probability':enc(risk),'uniform_hazard_integral':str(integrated_cap),
            'any_shortcut_jump_upper_bound':enc(cap),
            'winding_risk_corollary':'Any newly possible auxiliary winding requires a shortcut jump and is bounded by this event probability.'})
    result={'schema':'cyclecraft.v2.arb-finite-horizon.v1','status':'PASS'if all(c['pass']for c in checks)else'FAIL',
            'precision_bits':ctx.prec,'mpmath_reference_digits':mp.mp.dps,
            'arithmetic':'python-flint Arb rigorous ball matrix exponentials; exact rational source rates and dyadic endpoint serialization',
            'protocol':{'states':N,'edges':{k:v for k,v in GRAPHS.items()},'stage_duration':str(TAU),
                        'well_factors':[[str(v)for v in row]for row in W],
                        'core_conductances':[[str(v)for v in row]for row in QC],
                        'aux_conductance_cycle':'1/2','aux_conductance_star':'3/4','period':'6',
                        'observation_periods':10,'aux_counted_edge_index':3,'target_counted_edge_index':0},
            'serialization':'Exact dyadic lower/upper endpoint pairs are authoritative; float fields are nearest-binary64 display approximations only.',
            'certificates':out,'rare_shortcuts':shortcuts,'checks':checks,'elapsed_seconds':time.perf_counter()-st}
    (ROOT/'precision_cert_results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'status':result['status'],'checks':checks,'elapsed_seconds':result['elapsed_seconds'],
        'summary':{k:{'target_tail':v['target_lower_tail']['certified_probability_upper'],
                      'aux_survival_b1':v['auxiliary_survival']['1']}for k,v in out.items()}},indent=2))
    assert result['status']=='PASS'

if __name__=='__main__':main()
