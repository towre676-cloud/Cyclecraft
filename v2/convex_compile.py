"""Finite-speed occupation synthesis with fixed continuous gate conductances.

q and x are piecewise linear. Their product is quadratic and p is cubic.
The continuous-time dual certificate is implemented separately.
"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent

def graph(kind):
    aux={'star':[(0,3),(0,4)],'chain':[(0,3),(3,4)],'cycle':[(0,3),(3,4),(4,0)]}[kind]
    edges=[(0,1),(1,2),(2,0)]+aux;D=np.zeros((5,len(edges)))
    for e,(i,j)in enumerate(edges):D[i,e]=-1;D[j,e]=1
    return edges,D

def q_schedule(kind,N,amplitude=.625):
    edges,D=graph(kind);q=np.ones((N,len(edges)));q[:,0]=np.round((1+amplitude*np.sin(2*np.pi*np.arange(N)/N))*1024)/1024
    q[:,3:]=1 if kind!='cycle'else 2/3
    # Cycle 2/3 is replaced by nearest dyadic; the installation sum is recorded.
    q[:,3:]=np.round(q[:,3:]*1024)/1024
    if kind=='cycle':q[:,3:]=np.array([683,683,682])/1024
    return q

def partial_matrices(D,q0,q1,dt,s):
    # Integral_0^s q(r)*[(1-r)x0+r*x1] dt.
    dq=q1-q0
    a=dt*(q0*(s-s*s/2)+dq*(s*s/2-s**3/3))
    b=dt*(q0*s*s/2+dq*s**3/3)
    return (D*a)@D.T,(D*b)@D.T,a,b

def response_matrices(D,q0,q1,dt,s,mode):
    active=list(range(3))if mode=='baseline'else list(range(5))
    L0=(D*q0)@D.T;Ld=(D*(q1-q0))@D.T
    R=np.zeros((5,5))
    if mode=='baseline':R[:,3:]=-L0[:,3:]
    z=np.zeros((20,20));z[:5,:5]=R
    for i in range(3):z[i*5:(i+1)*5,(i+1)*5:(i+2)*5]=np.eye(5)
    ee=expm(z*(dt*s));F0=ee[:5,5:10];F1=ee[:5,10:15]/dt;F2=2*ee[:5,15:20]/dt**2
    B0=-L0[:,active];Bd=-Ld[:,active]
    C0=F0@B0+F1@(Bd-B0)-F2@Bd;C1=F1@B0+F2@Bd
    return ee[:5,:5],C0,C1

def solve(kind='star',mode='assisted',N=32,T=8.,S=.15,A=16.,Ebox=1.2,rmax=6.,pmin=.015):
    import cvxpy as cp
    edges,D=graph(kind);q=q_schedule(kind,N);dt=T/N
    p=cp.Variable((N,5));nu=3 if mode=='baseline'else 5;x=cp.Variable((N,nu));cons=[];dyn=[]
    cons+=[cp.sum(p,axis=1)==1,p>=pmin,x>=1e-5]
    quadrature=leggauss(3);qs=(quadrature[0]+1)/2;qw=quadrature[1]/2
    entropy=0;activity=0;objective=0;meanp=0
    for k in range(N):
        kp=(k+1)%N;Am,Bm,a,b=partial_matrices(D,q[k],q[kp],dt,1)
        Ph,C0,C1=response_matrices(D,q[k],q[kp],dt,1,mode)
        c=p[kp]-Ph@p[k]-C0@x[k]-C1@x[kp]==0;cons.append(c);dyn.append(c)
        objective+=a[0]*(x[k,0]-x[k,1])+b[0]*(x[kp,0]-x[kp,1])
        for s,w in zip(qs,qw):
            qv=(1-s)*q[k]+s*q[kp];Ph,C0,C1=response_matrices(D,q[k],q[kp],dt,s,mode)
            pv=Ph@p[k]+C0@x[k]+C1@x[kp];uv=(1-s)*x[k]+s*x[kp]
            xv=cp.hstack([uv,pv[3:]])if mode=='baseline'else uv
            meanp+=dt*w*pv;activity+=dt*w*cp.sum(cp.multiply(np.abs(D)@qv,xv))
            cons +=[pv>=pmin,xv>=np.exp(-Ebox)*pv,xv<=np.exp(Ebox)*pv,
                    cp.multiply(np.abs(D)@qv,xv)<=rmax*pv]
            for e,(i,j)in enumerate(edges):
                entropy+=dt*w*qv[e]*(cp.rel_entr(xv[i],xv[j])+cp.rel_entr(xv[j],xv[i]))
        xk=cp.hstack([x[k],p[k,3:]])if mode=='baseline'else x[k]
        cons +=[xk>=np.exp(-Ebox)*p[k],xk<=np.exp(Ebox)*p[k],cp.multiply(np.abs(D)@q[k],xk)<=rmax*p[k]]
    ce=entropy<=S;ca=activity<=A;cm=meanp==np.ones(5)*T/5;cons +=[ce,ca,cm]
    problem=cp.Problem(cp.Maximize(objective),cons)
    st=time.perf_counter();problem.solve(solver='CLARABEL',tol_gap_abs=2e-9,tol_feas=2e-9,tol_gap_rel=2e-9,max_iter=250)
    out={'architecture':kind,'mode':mode,'N':N,'period':T,'entropy_budget':S,'activity_budget':A,
         'energy_bound':Ebox,'escape_bound':rmax,'probability_floor':pmin,'status':problem.status,
         'objective':float(problem.value),'entropy_discrete':float(entropy.value),'activity':float(activity.value),
         'p_knots':p.value.tolist()if p.value is not None else None,'x_knots':x.value.tolist()if x.value is not None else None,
         'q_knots':q.tolist(),'dynamic_duals':[c.dual_value.tolist()for c in dyn],
         'entropy_dual':float(ce.dual_value),'activity_dual':float(ca.dual_value),
         'mean_occupation_dual':cm.dual_value.tolist(),'mean_occupation':(meanp.value/T).tolist(),
         'solver_seconds':time.perf_counter()-st,'edges':edges}
    return out

def evaluate(sol,phases):
    kind=sol['architecture'];edges,D=graph(kind);q=np.asarray(sol['q_knots']);x=np.asarray(sol['x_knots']);p=np.asarray(sol['p_knots']);N=len(q);dt=sol['period']/N
    arr=[]
    for phase in phases:
        z=(phase%1)*N;k=min(int(z),N-1);s=z-k;kp=(k+1)%N
        A,B,_,_=partial_matrices(D,q[k],q[kp],dt,s)
        if sol['mode']=='baseline':
            tau=dt*s;La=((D*q[k])@D.T)[3:,3:];ev,ee=np.linalg.eigh(La)
            bb=-((D*q[k])@D.T)[3:,0];inp=ee.T@bb;ini=ee.T@p[k,3:];udot=(x[kp,0]-x[k,0])/dt
            ex=np.exp(-ev*tau);one=-np.expm1(-ev*tau)
            aux=ee@(ex*ini+inp*(x[k,0]*one/ev+udot*(tau/ev-one/ev**2)))
            pv=np.zeros(5);pv[1:3]=(p[k]-A[:,:3]@x[k]-B[:,:3]@x[kp])[1:3];pv[3:]=aux;pv[0]=1-pv[1:].sum()
            xv=np.r_[(1-s)*x[k]+s*x[kp],aux]
        else:pv=p[k]-A@x[k]-B@x[kp];xv=(1-s)*x[k]+s*x[kp]
        qv=(1-s)*q[k]+s*q[kp]
        E=np.log(xv/pv);j=-qv*(D.T@xv);rates=qv[:,None]*(-D.T)*np.exp(E)[None,:]
        K=D@rates
        ell=np.array([(xv[i]-xv[jj])*np.log(xv[i]/xv[jj]) for i,jj in edges])
        entropy=np.sum(qv*ell);act=np.sum((np.abs(D)@qv)*xv)
        # All values are continuous; derivative may jump at knots.
        dx=(x[kp]-x[k])/dt;dp=-D@(qv*(D.T@xv))
        if sol['mode']=='baseline':dx=np.r_[dx,dp[3:]]
        slew=dx/xv-dp/pv
        arr.append((pv,xv,qv,E,j,entropy,act,K,slew))
    return arr

def validate(sol,samples=2049):
    edges,D=graph(sol['architecture'])
    phases=np.linspace(0,1,samples);rows=evaluate(sol,phases);p=np.array([r[0]for r in rows]);x=np.array([r[1]for r in rows]);q=np.array([r[2]for r in rows]);E=np.array([r[3]for r in rows]);j=np.array([r[4]for r in rows]);en=np.array([r[5]for r in rows]);ac=np.array([r[6]for r in rows]);slew=np.array([r[8]for r in rows]);T=sol['period']
    # Reintegrate the original master equation under the reconstructed controls.
    def rhs(t,y):return evaluate(sol,[t/T])[0][7]@y
    ss=solve_ivp(rhs,(0,T),p[0],method='DOP853',rtol=3e-12,atol=2e-14,t_eval=phases*T,max_step=T/len(sol['q_knots'])/2)
    counts=np.zeros_like(j)
    for k in range(1,len(phases)):counts[k]=counts[k-1]+(j[k]+j[k-1])*(T/(samples-1)/2)
    kk=len(sol['q_knots'])-1;pp=np.asarray(sol['p_knots']);xx=np.asarray(sol['x_knots']);qq=np.asarray(sol['q_knots']);Ph,C0,C1=response_matrices(D,qq[kk],qq[0],T/len(qq),1,sol['mode'])
    pend=Ph@pp[kk]+C0@xx[kk]+C1@xx[0]
    out={'min_probability':float(p.min()),'max_abs_energy':float(abs(E).max()),'max_escape':float(max(np.diag(-r[7]).max()for r in rows)),
         'entropy_trapezoid':float(np.trapezoid(en,phases*T)),'activity_trapezoid':float(np.trapezoid(ac,phases*T)),
         'max_abs_energy_slew':float(abs(slew).max()),'max_occupation_reintegration_error':float(abs(ss.y.T-p).max()),
         'knot_periodic_residual':float(max(abs(np.array(sol['p_knots'][0])-pend))),
         'aux_fixed_well_error':float(abs(E[:,3:]).max())if sol['mode']=='baseline'else None,
         'trajectory':{'phase':phases[::4].tolist(),'p':p[::4].tolist(),'x':x[::4].tolist(),'q':q[::4].tolist(),'E':E[::4].tolist(),'current':j[::4].tolist(),'counts':counts[::4].tolist()}}
    return out

if __name__=='__main__':
    N=int(sys.argv[1])if len(sys.argv)>1 else 32
    outputs=[]
    for kind in ['star','chain','cycle']:
        for mode in ['baseline','assisted']:
            r=solve(kind,mode,N=N);print(kind,mode,r['status'],r['objective'],r['entropy_discrete'],r['activity'],flush=True)
            r['validation']=validate(r)
            outputs.append(r)
            (ROOT/f'convex_raw_{N}.json').write_text(json.dumps(outputs,indent=2))
    (ROOT/f'convex_raw_{N}.json').write_text(json.dumps(outputs,indent=2))
