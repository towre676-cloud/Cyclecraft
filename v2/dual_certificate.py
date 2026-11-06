"""Continuous-time upper certificate for fixed-storage occupation control.

All proposal matrices and multipliers are binary64 dyadics. Arb encloses
continuous coefficients; LP dual residuals receive a rigorous compact-box
correction. No grid primal optimum is used as a continuous upper bound.
"""
import json,time,sys
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from flint import arb,ctx
from convex_compile import graph,evaluate
ctx.prec=100
ROOT=Path(__file__).resolve().parent
def ar(x):return arb(float(x))
def down(x):return np.nextafter(float(x.lower()),-np.inf)
def up(x):return np.nextafter(float(x.upper()),np.inf)
def rng(a,b):
    a=a if isinstance(a,arb)else ar(a);b=b if isinstance(b,arb)else ar(b)
    return arb((a+b)/2,abs(b-a)/2)

def lp_upper(c,A,b,AE,be,zmax):
    # Proposal only. Enclose its arbitrary dual and correct all residuals.
    rr=linprog(-c,A_ub=A,b_ub=b,A_eq=AE,b_eq=be,bounds=(0,None),method='highs')
    if not rr.success:raise RuntimeError(rr.message)
    y=np.maximum(-rr.ineqlin.marginals,0);z=-rr.eqlin.marginals
    ub=arb(0)
    for bi,yi in zip(b,y):ub+=ar(bi)*ar(yi)
    for bi,zi in zip(be,z):ub+=ar(bi)*ar(zi)
    residual=0.
    for j in range(len(c)):
        r=ar(c[j])
        for i in range(len(y)):r-=ar(A[i,j])*ar(y[i])
        for i in range(len(z)):r-=ar(AE[i,j])*ar(z[i])
        rp=max(up(r),0);ub+=ar(rp)*ar(zmax[j]);residual=max(residual,rp)
    return up(ub),residual,float(-rr.fun)

def certificate(sol,subdiv=64,phi_sign=-1,phi_style='endpoint'):
    assert sol['mode']=='baseline'
    edges,D=graph(sol['architecture']);q=np.array(sol['q_knots']);x=np.array(sol['x_knots']);p=np.array(sol['p_knots']);N=len(q);dt=sol['period']/N
    assert N>0 and N&(N-1)==0,'Certificate requires a power-of-two control partition'
    assert isinstance(subdiv,int)and subdiv>0 and subdiv&(subdiv-1)==0,'Certificate subdivisions must be a positive power of two'
    assert np.isfinite(dt)and dt>0,'Certificate period must be positive and finite'
    assert Fraction(dt)*N==Fraction(sol['period']),'Control interval must be an exact dyadic division of the source period'
    assert Fraction(dt/subdiv)*subdiv==Fraction(dt),'Certificate time weights must be exactly representable; avoid underflow'
    yd=np.array(sol['dynamic_duals'])
    phi=phi_sign*(np.roll(yd,1,axis=0)if phi_style=='endpoint'else(yd+np.roll(yd,1,axis=0))/2)
    beta=max(sol['entropy_dual'],0);gamma=max(sol['activity_dual'],0)
    lam=-np.array(sol['mean_occupation_dual'])
    eb=np.nextafter(sol['energy_bound'],np.inf)
    lower=down(ar(-eb).exp());upper=up(ar(eb).exp())
    xmin=ar(lower)*ar(sol['probability_floor']);xmax=ar(upper)
    jeffreys_cap=up((xmax-xmin)*(xmax/xmin).log())
    assert jeffreys_cap<64,'Increase epigraph cap for this energy/probability box'
    nv=10+len(edges)
    zmax=np.r_[np.ones(5),np.ones(5)*upper,np.ones(len(edges))*64]
    # Common polytope: nonnegative p,x; normalization; fixed auxiliary wells.
    AE=np.zeros((3,nv));AE[0,:5]=1;AE[1,3]=-1;AE[1,8]=1;AE[2,4]=-1;AE[2,9]=1;be=np.array([1.,0,0])
    staticA=[];staticb=[]
    for i in range(5):
        r=np.zeros(nv);r[i]=-1;staticA.append(r);staticb.append(-sol['probability_floor'])
        r=np.zeros(nv);r[i]=lower;r[5+i]=-1;staticA.append(r);staticb.append(0.)
        r=np.zeros(nv);r[i]=-upper;r[5+i]=1;staticA.append(r);staticb.append(0.)
    for e in range(len(edges)):
        r=np.zeros(nv);r[10+e]=1;staticA.append(r);staticb.append(64.)
    total=ar(beta)*ar(np.nextafter(sol['entropy_budget'],np.inf))+ar(gamma)*ar(sol['activity_budget'])
    for l in lam:total-=ar(l)*ar(sol['period'])/5
    maxres=0.;rawtotal=0.;st=time.perf_counter()
    for k in range(N):
        kp=(k+1)%N;dphi=(phi[kp]-phi[k])/dt
        for j in range(subdiv):
            a=j/subdiv;b=(j+1)/subdiv;mid=(a+b)/2
            qi=[rng(ar(q[k,e])+ar(a)*(ar(q[kp,e])-ar(q[k,e])),ar(q[k,e])+ar(b)*(ar(q[kp,e])-ar(q[k,e])))for e in range(len(edges))]
            ph=[rng(ar(phi[k,i])+ar(a)*(ar(phi[kp,i])-ar(phi[k,i])),ar(phi[k,i])+ar(b)*(ar(phi[kp,i])-ar(phi[k,i])))for i in range(5)]
            # Tangent points from exact baseline occupation response.
            row=evaluate(sol,[(k+mid)/N])[0];xx=row[1]
            deg=[arb(0)for _ in range(5)];cx=[arb(0)for _ in range(5)];planes=[];cz=[]
            for e,(u,v)in enumerate(edges):
                qa=qi[e];deg[u]+=qa;deg[v]+=qa
                tmp=qa*(ph[u]-ph[v]);cx[u]+=tmp;cx[v]-=tmp
                qm=ar(down(qa));cz.append(-ar(beta)*qm)
                for offset in [-4,-2,-1,-.5,-.25,-.125,0,.125,.25,.5,1,2,4]:
                    ratio=float(xx[u]/xx[v]*np.exp(offset));rr=ar(ratio);lg=rr.log();g0=lg+1-1/rr;g1=-lg+1-rr
                    plane=np.zeros(nv);plane[5+u]=down(g0);plane[5+v]=down(g1);plane[10+e]=-1;planes.append(plane)
            cx[0]+=qi[0];cx[1]-=qi[0]
            for i in range(5):cx[i]-=ar(gamma)*deg[i]
            cpcoef=[up(-(ar(phi[kp,i])-ar(phi[k,i]))/ar(dt)+ar(lam[i]))for i in range(5)]
            cc=np.r_[cpcoef,np.array([up(z)for z in cx]),np.array([up(z)for z in cz])]
            AA=list(staticA);bb=list(staticb)
            AA.extend(planes);bb.extend([0.]*len(planes))
            for i in range(5):
                rr=np.zeros(nv);rr[i]=-sol['escape_bound'];rr[5+i]=down(deg[i]);AA.append(rr);bb.append(0.)
            ub,res,raw=lp_upper(cc,np.array(AA),np.array(bb),AE,be,zmax)
            total+=ar(dt/subdiv)*ar(ub);rawtotal+=dt/subdiv*raw;maxres=max(maxres,res)
    return {'architecture':sol['architecture'],'continuous_upper':up(total),'subintervals':N*subdiv,'subdivisions_per_control_interval':subdiv,
            'entropy_multiplier':beta,'activity_multiplier':gamma,'mean_multiplier':lam.tolist(),'phi_sign':phi_sign,
            'phi_style':phi_style,'epigraph_upper':64.,'jeffreys_box_upper':jeffreys_cap,'epigraph_bound_reason':'x_i between exp(-energy_bound)*p_floor and exp(energy_bound)',
            'max_LP_dual_residual_correction_coefficient':maxres,'arb_precision_bits':ctx.prec,
            'method':'continuous occupation dual; interval coefficient bounds; Jeffreys entropy supporting planes; LP dual residual box correction',
            'seconds':time.perf_counter()-st}

if __name__=='__main__':
    N=int(sys.argv[1])if len(sys.argv)>1 else 32;sub=int(sys.argv[2])if len(sys.argv)>2 else 64
    sols=json.loads((ROOT/f'convex_raw_{N}.json').read_text());out=[]
    for s in sols:
        if s['mode']=='baseline':
            c=certificate(s,sub);out.append(c);print(json.dumps(c),flush=True)
    (ROOT/f'dual_bounds_{N}_{sub}.json').write_text(json.dumps(out,indent=2))
