"""Cyclecraft II: smooth periodic statistics, graph noise and count risk.

Column generators: dp/dt=Kp. D has -1 tail,+1 head. One orientation per
undirected edge. Monitored observables are signed crossings on edge 0
(productive triangle) and edge 3 (auxiliary channel). Smooth protocols are
harmonics in a local optimal friction-normalized plane. All auxiliary
architectures have two units of installed conductance; core q_e=1.

Run OPENBLAS_NUM_THREADS=1 python statistics.py.
"""
from pathlib import Path
from dataclasses import dataclass
import argparse,json,time,sys,platform
from types import SimpleNamespace
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq,minimize_scalar

ROOT=Path(__file__).resolve().parent
RTOL=2e-12
ATOL=3e-14
CORE=[(0,1),(1,2),(2,0)]
AUX={'cycle':[(0,3),(3,4),(4,0)],'star':[(0,3),(0,4)],'chain':[(0,3),(3,4)]}

def integrate(fun,span,y0,protocol,**kwargs):
    """Respect every known polynomial knot; harmonic paths use ordinary DOP853."""
    knots=getattr(protocol,'breakpoints',None)
    if knots is None:return solve_ivp(fun,span,y0,**kwargs)
    a,b=span;low,high=sorted(span);inside=[float(t)for t in knots if low<t<high]
    points=np.array([low]+sorted(inside)+[high]);forward=b>a
    intervals=list(zip(points[:-1],points[1:]));intervals=intervals if forward else [(j,i)for i,j in intervals[::-1]]
    opts=dict(kwargs);opts['dense_output']=True;sols=[];y=np.asarray(y0);nfev=0
    for aa,bb in intervals:
        lo,hi=sorted((aa,bb))
        # One-sided values at a knot keep jump-discontinuous work/slew terms
        # on the current polynomial span. The generator itself is continuous.
        def local_fun(t,z,lo=lo,hi=hi):
            if t<=lo:t=np.nextafter(lo,hi)
            if t>=hi:t=np.nextafter(hi,lo)
            return fun(t,z)
        sol=solve_ivp(local_fun,(aa,bb),y,**opts)
        if not sol.success:return sol
        sols.append(sol);y=sol.y[:,-1];nfev+=sol.nfev
    def dense(t):
        def one(tt):
            idx=int(np.clip(np.searchsorted(points,tt,side='right')-1,0,len(sols)-1))
            return sols[idx if forward else len(sols)-1-idx].sol(tt)
        tt=np.asarray(t)
        if tt.ndim==0:return one(float(tt))
        return np.column_stack([one(float(v))for v in tt])
    return SimpleNamespace(success=True,y=np.column_stack([y0,y]),nfev=nfev,sol=dense)

@dataclass
class Model:
    architecture:str
    def __post_init__(self):
        self.n=5;self.edges=CORE+AUX[self.architecture];self.m=len(self.edges)
        self.D=np.zeros((5,self.m))
        for e,(i,j) in enumerate(self.edges):self.D[i,e]=-1;self.D[j,e]=1
        self.q0=np.r_[np.ones(3),np.ones(self.m-3)*2/(self.m-3)]
        self.U=np.vstack([np.zeros(4),np.eye(4)])
        self.V=np.zeros((self.m,4));self.V[0,1]=1
        self.count=np.zeros((self.m,2));self.count[0,0]=1;self.count[3,1]=1
    def matrices(self,lam):
        E=self.U@lam;q=self.q0*np.exp(-self.V@lam);w=np.exp(E)
        C=-q[:,None]*self.D.T*w[None,:];K=self.D@C
        pi=np.exp(-E);pi/=pi.sum()
        return K,C,pi,E,q
    def geometry(self,lam=None):
        if lam is None:lam=np.zeros(4)
        K,C,pi,E,q=self.matrices(lam);L=(self.D*q)@self.D.T
        Lp=np.linalg.inv(L+np.ones((5,5))/5)-np.ones((5,5))/5
        H=q[:,None]*self.D.T@Lp;P=np.eye(self.m)-H@self.D
        dpi=-pi[:,None]*(self.U-pi@self.U)
        G=np.exp(-E).sum()*dpi.T@Lp@dpi
        F=np.zeros((self.m,4,4))
        for a in range(4):
            dH=-P@(self.V[:,a,None]*H)
            for b in range(4):F[:,a,b]+=dH@dpi[:,b];F[:,b,a]-=dH@dpi[:,b]
        return G,F,H,Lp

def top_plane(F,G):
    val,vec=eigh(G);S=(vec/np.sqrt(val))@vec.T;A=S@F@S
    sigma=np.sqrt(max(eigh(-A@A)[0][-1],0));Proj=-A@A/sigma**2
    a=Proj[:,0]
    if np.linalg.norm(a)<1e-10:a=Proj[:,np.argmax(np.diag(Proj))]
    a/=np.linalg.norm(a);b=-A@a/sigma
    return S@a,S@b,float(sigma)

@dataclass
class Harmonic:
    model:Model
    u:np.ndarray
    v:np.ndarray
    amplitude:float
    period:float
    def __call__(self,t):
        th=2*np.pi*t/self.period
        lam=self.amplitude*(self.u*np.cos(th)+self.v*np.sin(th))
        return self.model.matrices(lam)
    def velocity(self,t):
        th=2*np.pi*t/self.period
        return 2*np.pi*self.amplitude/self.period*(-self.u*np.sin(th)+self.v*np.cos(th))
    def energy_velocity(self,t):
        return self.model.U@self.velocity(t)
    def coordinates(self,t):
        th=2*np.pi*t/self.period
        return self.amplitude*(self.u*np.cos(th)+self.v*np.sin(th))
    def rate_dominator(self):
        ranges=self.amplitude*np.sqrt(self.u**2+self.v**2)
        emax=self.model.U@ranges
        qmax=self.model.q0*np.exp(np.abs(self.model.V)@ranges)
        exits=np.zeros(5)
        for e,(i,j) in enumerate(self.model.edges):exits[i]+=qmax[e]*np.exp(emax[i]);exits[j]+=qmax[e]*np.exp(emax[j])
        return float(exits.max()*(1+1e-13))

def stationary(M):
    A=np.eye(5)-M;A[-1]=1;b=np.zeros(5);b[-1]=1
    return np.linalg.solve(A,b)

def periodic_solution(protocol,samples=301,rtol=RTOL,atol=ATOL):
    n=5;mod=protocol.model;m=mod.m;T=protocol.period
    sm=integrate(lambda t,y:(protocol(t)[0]@y.reshape(n,n)).ravel(),(0,T),np.eye(n).ravel(),protocol,method='DOP853',rtol=rtol,atol=atol)
    assert sm.success
    M=sm.y[:,-1].reshape(n,n);p0=stationary(M)
    def fun(t,y):
        K,C,pi,E,q=protocol(t);p=y[:n];j=C@p;x=np.exp(E)*p
        lgr=np.log(x[np.array(mod.edges)[:,0]]/x[np.array(mod.edges)[:,1]])
        sigma=j@lgr;work=p@protocol.energy_velocity(t)
        traffic=np.array([q[e]*(x[i]+x[k])for e,(i,k)in enumerate(mod.edges)])
        return np.r_[K@p,j,sigma,work,traffic.sum()]
    ss=integrate(fun,(0,T),np.r_[p0,np.zeros(m+3)],protocol,method='DOP853',rtol=rtol,atol=atol,dense_output=True)
    assert ss.success
    ts=np.linspace(0,T,samples);ys=ss.sol(ts);yend=ss.y[:,-1]
    pp=ys[:n].T;curr=[];pis=[];Es=[];qs=[];lams=[];maxexit=0.;holoerr=0.;floor=[]
    for t,p in zip(ts,pp):
        K,C,pi,E,q=protocol(t);curr.append((C@p).tolist());pis.append(pi.tolist());Es.append(E.tolist());qs.append(q.tolist());lams.append(protocol.coordinates(t).tolist());maxexit=max(maxexit,float((-np.diag(K)).max()))
        L=(mod.D*q)@mod.D.T;Lp=np.linalg.inv(L+np.ones((5,5))/5)-np.ones((5,5))/5;H=q[:,None]*mod.D.T@Lp
        holoerr=max(holoerr,float(abs(H@(K@p)-C@p).max()))
    out={'amplitude':float(protocol.amplitude),'period':T,'mean_target':float(yend[n]),'mean_aux':float(yend[n+3]),'edge_counts':yend[n:n+m].tolist(),
         'entropy':float(yend[n+m]),'work':float(yend[n+m+1]),'activity':float(yend[n+m+2]),'max_exit_rate_sampled':maxexit,
         'analytic_rate_dominator':protocol.rate_dominator(),'work_entropy_error':float(abs(yend[n+m]-yend[n+m+1])),
         'periodic_residual':float(abs(yend[:n]-p0).max()),'mass_residual':float(abs(pp.sum(axis=1)-1).max()),
         'current_divergence':float(abs(mod.D@yend[n:n+m]).max()),'exact_holonomy_residual':holoerr,
         'trajectory':{'phase':(ts/T).tolist(),'p':pp.tolist(),'pi':pis,'current':curr,'counts':ys[n:n+m].T.tolist(),'E':Es,'q':qs,'lambda':lams,
                       'entropy':ys[n+m].tolist(),'activity':ys[n+m+2].tolist()},'initial_probability':p0.tolist(),'nfev':sm.nfev+ss.nfev}
    return out,M,p0,ss

def tilt_derivatives(protocol):
    T=protocol.period;mod=protocol.model
    # P,R_L,R_R,S_LL,S_LR,S_RR are derivatives of the period propagator.
    def fun(t,y):
        K=protocol(t)[0];P,RL,RR,SLL,SLR,SRR=y.reshape(6,5,5)
        K1=[];K2=[]
        for e in [0,3]:
            i,j=mod.edges[e];A=np.zeros((5,5));B=A.copy()
            A[j,i]=K[j,i];A[i,j]=-K[i,j];B[j,i]=K[j,i];B[i,j]=K[i,j]
            K1.append(A);K2.append(B)
        A,B=K1;A2,B2=K2
        return np.array([K@P,K@RL+A@P,K@RR+B@P,K@SLL+2*A@RL+A2@P,
                         K@SLR+A@RR+B@RL,K@SRR+2*B@RR+B2@P]).ravel()
    ini=np.zeros((6,5,5));ini[0]=np.eye(5)
    s=integrate(fun,(0,T),ini.ravel(),protocol,method='DOP853',rtol=RTOL,atol=ATOL)
    assert s.success
    return s.y[:,-1].reshape(6,5,5),s.nfev

def counting_statistics(blocks,p0,horizons=(1,10,100)):
    M,RL,RR,SLL,SLR,SRR=blocks;ones=np.ones(5);Rs=[RL,RR];Ss=[[SLL,SLR],[SLR,SRR]]
    mu=np.array([ones@A@p0 for A in Rs]);correctors=[]
    for a in range(2):correctors.append(np.linalg.solve(np.eye(5)-M+np.outer(p0,ones),(Rs[a]-mu[a]*np.eye(5))@p0))
    cov=np.zeros((2,2))
    for a in range(2):
        for b in range(2):cov[a,b]=ones@(Ss[a][b]@p0+Rs[a]@correctors[b]+Rs[b]@correctors[a])-mu[a]*mu[b]
    state=np.zeros((6,5));state[0]=p0;finite=[]
    for n in range(1,max(horizons)+1):
        P,a,b,aa,ab,bb=state
        state=np.array([M@P,RL@P+M@a,RR@P+M@b,SLL@P+2*RL@a+M@aa,SLR@P+RL@b+RR@a+M@ab,SRR@P+2*RR@b+M@bb])
        if n in horizons:
            means=np.array([ones@state[1],ones@state[2]])
            raw2=np.array([[ones@state[3],ones@state[4]],[ones@state[4],ones@state[5]]]);cv=raw2-np.outer(means,means)
            finite.append({'cycles':n,'mean_target':float(means[0]),'mean_aux':float(means[1]),'variance_target':float(cv[0,0]),'variance_aux':float(cv[1,1]),'covariance':cv.tolist()})
    return {'mean':mu.tolist(),'covariance_asymptotic':cov.tolist(),'variance_target':float(cov[0,0]),'variance_aux':float(cov[1,1]),
            'cross_covariance':float(cov[0,1]),'asymptotic_eigenvalues':eigh(cov)[0].tolist(),
            'copies_cycles_for_snr_5':float(25*cov[0,0]/mu[0]**2),'finite':finite}

def poisson_check(protocol,M,p0,probsol,mu):
    T=protocol.period;mod=protocol.model;one=np.ones(5);drift=np.array(mu)/T
    def rhs(t,y):
        K=protocol(t)[0];h=y.reshape(5,2);f=np.zeros((5,2))
        for a,e in enumerate([0,3]):
            i,j=mod.edges[e];f[i,a]+=K[j,i];f[j,a]-=K[i,j]
        return (-K.T@h+one[:,None]*drift-f).ravel()
    s0=integrate(rhs,(T,0),np.zeros(10),protocol,method='DOP853',rtol=RTOL,atol=ATOL)
    b=s0.y[:,-1].reshape(5,2);b-=one[:,None]*(p0@b)[None,:]
    h0=np.linalg.solve(np.eye(5)-M.T+np.outer(one,p0),b)
    sh=integrate(rhs,(T,0),h0.ravel(),protocol,method='DOP853',rtol=RTOL,atol=ATOL,dense_output=True)
    def quadrature(t,y):
        K,C,pi,E,q=protocol(t);p=probsol.sol(t)[:5];x=np.exp(E)*p
        traffic=np.array([q[e]*(x[i]+x[j])for e,(i,j)in enumerate(mod.edges)])
        h=sh.sol(t).reshape(5,2);inc=mod.count+mod.D.T@h
        A=np.diag(traffic);La=(mod.D*traffic)@mod.D.T
        Lap=np.linalg.inv(La+np.ones((5,5))/5)-np.ones((5,5))/5
        F=A-A@mod.D.T@Lap@mod.D@A
        cv=inc.T@A@inc;floor=mod.count.T@F@mod.count
        return np.r_[cv[0,0],cv[0,1],cv[1,1],floor[0,0],floor[0,1],floor[1,1]]
    sq=integrate(quadrature,(0,T),np.zeros(6),protocol,method='DOP853',rtol=RTOL,atol=ATOL)
    a,b,c,d,e,f=sq.y[:,-1];cv=np.array([[a,b],[b,c]]);floor=np.array([[d,e],[e,f]])
    return {'covariance':cv.tolist(),'graph_floor':floor.tolist(),'excess_eigenvalues':eigh(cv-floor)[0].tolist(),
            'corrector_periodic_residual':float(abs(sh.y[:,-1].reshape(5,2)-h0).max())}

def safe_monodromy(protocol,bound):
    """Exact killed count lift for signed auxiliary count |N_R|<=bound."""
    T=protocol.period;mod=protocol.model;size=5*(2*bound+1);monitor=3
    # Integrate sparse-action columns; transitions leaving count range are killed.
    def fun(t,y):
        K=protocol(t)[0];Y=y.reshape(2*bound+1,5,size);out=np.einsum('ij,kjl->kil',K,Y)
        i,j=mod.edges[monitor];f=K[j,i];r=K[i,j]
        # Remove within-count monitor transitions and place across count slices.
        out[:,j]-=f*Y[:,i];out[:,i]-=r*Y[:,j]
        out[1:,j]+=f*Y[:-1,i];out[:-1,i]+=r*Y[1:,j]
        return out.ravel()
    sol=integrate(fun,(0,T),np.eye(size).ravel(),protocol,method='DOP853',rtol=5e-12,atol=5e-14)
    assert sol.success
    return sol.y[:,-1].reshape(size,size)

def survival(protocol,p0,bounds=(1,2),horizons=(1,10,100)):
    out=[]
    for b in bounds:
        if protocol.model.architecture!='cycle':
            # For a bridge N_R is exactly an endpoint indicator difference.
            for n in horizons:out.append({'bound':b,'cycles':n,'survival':1.,'failure':0.,'method':'exact bridge endpoint identity'})
            continue
        M=safe_monodromy(protocol,b);v=np.zeros(M.shape[0]);v[b*5:(b+1)*5]=p0
        for n in range(1,max(horizons)+1):
            v=M@v
            if n in horizons:
                s=float(v.sum());out.append({'bound':b,'cycles':n,'survival':s,'failure':float(1-s),'method':'killed count-lift propagator'})
    return out

def stochastic_path(protocol,p0,seed,duration_cycles=10):
    """Exact thinning using a proved analytic bound on all state exit rates."""
    rng=np.random.default_rng(seed);R=protocol.rate_dominator();duration=protocol.period*duration_cycles
    t=0.;i=int(rng.choice(5,p=p0));target=0;aux=0;ts=[0.];states=[i];ll=[0];rr=[0];proposals=0
    while True:
        t+=rng.exponential(1/R);proposals+=1
        if t>duration:break
        K=protocol(t)[0];rates=K[:,i].copy();rates[i]=0;ex=rates.sum();assert ex<=R*(1+1e-12)
        if rng.random()*R>=ex:continue
        j=int(rng.choice(5,p=rates/ex));
        for a,e in enumerate([0,3]):
            u,v=protocol.model.edges[e];inc=1 if (i,j)==(u,v) else -1 if(i,j)==(v,u)else 0
            if a==0:target+=inc
            else:aux+=inc
        i=j;ts.append(float(t));states.append(i);ll.append(target);rr.append(aux)
    ts.append(duration);states.append(i);ll.append(target);rr.append(aux)
    return {'seed':seed,'duration':duration,'cycles':duration_cycles,'time':ts,'state':states,'target_count':ll,'aux_count':rr,'rate_dominator':R,'proposals':proposals,
            'method':'exact time-dependent thinning with analytic harmonic-rate dominator'}

def mgf_monodromy(protocol,theta):
    mod=protocol.model;T=protocol.period
    def fun(t,y):
        K=protocol(t)[0].copy()
        for a,e in enumerate([0,3]):
            i,j=mod.edges[e];K[j,i]*=np.exp(theta[a]);K[i,j]*=np.exp(-theta[a])
        return (K@y.reshape(5,5)).ravel()
    s=integrate(fun,(0,T),np.eye(5).ravel(),protocol,method='DOP853',rtol=RTOL,atol=ATOL)
    assert s.success
    return s.y[:,-1].reshape(5,5)

def finite_chernoff(protocol,p0,ncycles=10):
    """A finite-horizon lower-tail bound for target net count <=0."""
    cache={}
    def logbound(s):
        if s not in cache:
            M=mgf_monodromy(protocol,[-s,0]);v=np.linalg.matrix_power(M,ncycles)@p0
            cache[s]=float(np.log(v.sum()))
        return cache[s]
    opt=minimize_scalar(logbound,bounds=(1e-5,.15),method='bounded',options={'xatol':1e-5})
    return {'cycles':ncycles,'event':'target signed count <= 0','theta':float(-opt.x),'log_bound_per_copy':float(opt.fun),'bound':float(min(1,np.exp(opt.fun))),
            'independent_ensembles':[{'copies':m,'bound':float(min(1,np.exp(m*opt.fun)))}for m in [1,100,1000,10000]],
            'method':'finite-period monodromy Chernoff bound, numerical binary64'}

def analyze_case(protocol,design,full=True,seed=20261031):
    r,M,p0,ps=periodic_solution(protocol)
    r.update({'architecture':protocol.model.architecture,'design':design,'edges':[list(e)for e in protocol.model.edges]})
    if not full:return r
    blocks,nfev=tilt_derivatives(protocol);cs=counting_statistics(blocks,p0);r.update(cs)
    pc=poisson_check(protocol,M,p0,ps,cs['mean']);r['poisson_check']=pc
    r['floquet_poisson_error']=float(abs(np.array(pc['covariance'])-np.array(cs['covariance_asymptotic'])).max())
    r['counting_mean_error']=float(max(abs(cs['mean'][0]-r['mean_target']),abs(cs['mean'][1]-r['mean_aux'])))
    r['auxiliary_is_bridge']=protocol.model.architecture!='cycle'
    r['asymptotic_count_covariance_rank_exact']=1 if r['auxiliary_is_bridge']else 2
    if r['auxiliary_is_bridge']:
        # Independent finite-horizon variance from the endpoint indicator law.
        descendants={3}if protocol.model.architecture=='star'else {3,4}
        mask=np.array([i in descendants for i in range(5)],dtype=float);prob=float(mask@p0);errs=[]
        for row in r['finite']:
            Pn=np.linalg.matrix_power(M,row['cycles'])
            endpoint_var=2*prob-2*float(mask@Pn@(mask*p0))
            row['endpoint_variance_aux']=endpoint_var
            errs.append(abs(endpoint_var-row['variance_aux']))
        r['bridge_endpoint_variance_error']=float(max(errs))
    r['survival']=survival(protocol,p0)
    for row in r['finite']:
        row['survival_aux']=next(x['survival']for x in r['survival']if x['bound']==1 and x['cycles']==row['cycles'])
    r['sample_path']=stochastic_path(protocol,p0,seed)
    r['chernoff']=finite_chernoff(protocol,p0)
    r['nfev']+=nfev
    return r

class RepairedProtocol:
    """Actual rate schedule reconstructed from a certified rational p,x,q path."""
    def __init__(self,certificate):
        from primal_certificate import evaluate_protocol
        self.evaluate=evaluate_protocol;self.certificate=certificate
        self.source=certificate['protocol'];self.model=Model(self.source['architecture'])
        self.period=float(self.source['period']);self.amplitude=float(63/64)
        self.breakpoints=np.linspace(0,self.period,self.source['N']+1)
    def state(self,t):return self.evaluate(self.source,(t/self.period)%1)
    def __call__(self,t):
        a=self.state(t);return a['K'],a['C'],a['pi'],a['E'],a['q']
    def energy_velocity(self,t):return self.state(t)['E_dot']
    def coordinates(self,t):return self.state(t)['E'][1:]
    def rate_dominator(self):
        # This cap was validated for the exact source; floating replay has slack.
        return float(self.certificate['budgets']['escape_rate']*(1+1e-12))

def compiled_main(path):
    start=time.perf_counter();raw=json.loads(Path(path).read_text());cases=[];checks=[]
    for index,cert in enumerate(raw['results']):
        if cert['status']!='PASS':raise ValueError('Refusing an uncertified source protocol')
        print('compiled architecture',cert['architecture'],flush=True)
        protocol=RepairedProtocol(cert);r=analyze_case(protocol,'certified_assisted',seed=20261120+index)
        r['control_parameterization']='certified rational affine x and q; cubic p; E=log(x/p)'
        r['construction_shrink_factor']=cert['shrink_factor'];r.pop('amplitude',None)
        r['source_certificate_status']=cert['status'];r['source_N']=cert['input_N']
        r['certified_mean_target']=cert['objective'];r['certified_activity']=cert['activity'];r['entropy_enclosure']=cert['entropy_enclosure']
        r['source_protocol_file']=Path(path).name
        wanted=np.array([protocol.state(phase*protocol.period)['p']for phase in r['trajectory']['phase']])
        r['occupation_replay_error']=float(abs(np.asarray(r['trajectory']['p'])-wanted).max())
        r['certified_mean_replay_error']=float(abs(r['mean_target']-cert['objective']))
        r['certified_activity_replay_error']=float(abs(r['activity']-cert['activity']))
        tests=[('occupation replay',r['occupation_replay_error'],1e-10),('certified mean',r['certified_mean_replay_error'],1e-10),
               ('certified activity',r['certified_activity_replay_error'],1e-10),('work entropy',r['work_entropy_error'],1e-10),
               ('Floquet Poisson',r['floquet_poisson_error'],2e-9),('periodicity',r['periodic_residual'],1e-10)]
        if r['auxiliary_is_bridge']:
            tests+=[('bridge diffusion',abs(r['variance_aux']),1e-9),('endpoint variance',r['bridge_endpoint_variance_error'],1e-9)]
        else:tests+=[('positive auxiliary noise',-r['variance_aux'],0.)]
        for label,value,bound in tests:checks.append({'name':cert['architecture']+' '+label,'value':value,'bound':bound,'pass':bool(value<bound)})
        cases.append(r)
        print(cert['architecture'],'mean',r['mean_target'],'Vtarget',r['variance_target'],'mean error',r['certified_mean_replay_error'],flush=True)
    out={'schema':'cyclecraft-v2-certified-protocol-statistics-v1','source_file':Path(path).name,'precision':'Binary64 independent piecewise DOP853 and analytic counting derivatives; source protocol feasibility and mean are separately exact/interval certified.',
         'cases':cases,'checks':checks,'all_checks_pass':all(x['pass']for x in checks),'elapsed_seconds':time.perf_counter()-start}
    (ROOT/'statistics_compiled_results.json').write_text(json.dumps(out,indent=2));assert out['all_checks_pass'];print('COMPILED FINAL',len(checks),'PASS',flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--quick',action='store_true');parser.add_argument('--compiled');args=parser.parse_args();start=time.perf_counter()
    if args.compiled:return compiled_main(args.compiled)
    output={'schema':'cyclecraft-v2-statistics-v1','precision':'IEEE754 binary64 DOP853; analytic tilted derivatives and independent periodic Poisson corrector',
            'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'models':{},'cases':[],'checks':[]}
    def check(name,value,bound):output['checks'].append({'name':name,'value':float(value),'bound':bound,'pass':bool(value<bound)})
    for ai,name in enumerate(['cycle','star','chain']):
        print('architecture',name,flush=True);mod=Model(name);G,F,H,Lp=mod.geometry();Ft=F[:3].mean(axis=0)
        ub,vb,kb=top_plane(Ft[:2,:2],G[:2,:2]);ub=np.r_[ub,0,0];vb=np.r_[vb,0,0]
        ua,va,ka=top_plane(Ft,G)
        modelout={'edges':[list(e)for e in mod.edges],'D':mod.D.tolist(),'base_q':mod.q0.tolist(),'installed_aux_conductance':float(mod.q0[3:].sum()),
                  'geometry':{'G':G.tolist(),'F_target':Ft.tolist(),'kappa_baseline':kb,'kappa_assisted':ka,'gain':ka/kb-1,
                              'baseline_u':ub.tolist(),'baseline_v':vb.tolist(),'assisted_u':ua.tolist(),'assisted_v':va.tolist()},'cases':{}}
        T=20.;amp=.12
        base=Harmonic(mod,ub,vb,amp,T);ass=Harmonic(mod,ua,va,amp,T)
        rb=analyze_case(base,'baseline',full=not args.quick,seed=20261031+ai*10)
        ra=analyze_case(ass,'assisted',full=not args.quick,seed=20261032+ai*10)
        entropy_target=rb['entropy']
        def discrepancy(a):return periodic_solution(Harmonic(mod,ua,va,a,T),samples=11,rtol=1e-11,atol=1e-13)[0]['entropy']-entropy_target
        match=brentq(discrepancy,.005,.3,xtol=2e-13)
        rm=analyze_case(Harmonic(mod,ua,va,match,T),'assisted_matched',full=not args.quick,seed=20261033+ai*10)
        modelout['matched_entropy_gain']=rm['mean_target']/rb['mean_target']-1
        modelout['matched_entropy_amplitude']=match
        for label,r in [('baseline',rb),('assisted',ra),('assisted_matched',rm)]:
            modelout['cases'][label]=r;output['cases'].append(r)
            check(name+' '+label+' mass',r['mass_residual'],1e-10)
            check(name+' '+label+' periodicity',r['periodic_residual'],1e-10)
            check(name+' '+label+' work entropy',r['work_entropy_error'],1e-10)
            check(name+' '+label+' aux mean',abs(r['mean_aux']),1e-10)
            if not args.quick:
                check(name+' '+label+' counting mean',r['counting_mean_error'],1e-10)
                check(name+' '+label+' Floquet Poisson',r['floquet_poisson_error'],2e-9)
                if name!='cycle':
                    check(name+' '+label+' bridge diffusion',abs(r['variance_aux']),1e-9)
                    check(name+' '+label+' path bridge bound',max(abs(np.array(r['sample_path']['aux_count'])))-1,1e-13)
                    check(name+' '+label+' endpoint variance',r['bridge_endpoint_variance_error'],1e-9)
                else:check(name+' '+label+' positive aux noise',-r['variance_aux'],0)
        check(name+' matched entropy',abs(rm['entropy']-entropy_target),1e-10)
        if not args.quick:
            modelout['period_sweep']=[]
            for period in [2.,5.,10.,20.,40.,80.]:
                for design,u,v in [('baseline',ub,vb),('assisted',ua,va)]:
                    if period==T:rr=rb if design=='baseline'else ra
                    else:rr=periodic_solution(Harmonic(mod,u,v,amp,period),samples=41)[0]
                    modelout['period_sweep'].append({k:rr[k]for k in ['period','amplitude','mean_target','mean_aux','entropy','activity','max_exit_rate_sampled','periodic_residual']}|{'design':design})
            slowbase=Harmonic(mod,ub,vb,amp,80.);slowr,slowM,slowp,slowps=periodic_solution(slowbase,samples=41)
            target=slowr['entropy']
            def slow_discrepancy(a):return periodic_solution(Harmonic(mod,ua,va,a,80.),samples=11,rtol=1e-11,atol=1e-13)[0]['entropy']-target
            slowamp=brentq(slow_discrepancy,.005,.3,xtol=2e-13)
            slows=[]
            for design,pp,rr in [('baseline',slowbase,slowr),('assisted_matched',Harmonic(mod,ua,va,slowamp,80.),None)]:
                if rr is None:rr,_,p0slow,_=periodic_solution(pp,samples=41)
                else:p0slow=slowp
                blocks,_=tilt_derivatives(pp);cc=counting_statistics(blocks,p0slow)
                compact={k:rr[k]for k in ['period','amplitude','mean_target','mean_aux','entropy','activity','max_exit_rate_sampled','periodic_residual']}
                compact.update({'design':design,'variance_target':cc['variance_target'],'variance_aux':cc['variance_aux'],'covariance_asymptotic':cc['covariance_asymptotic'],'finite':cc['finite']})
                slows.append(compact)
            modelout['slow_matched']=slows;modelout['slow_matched_gain']=slows[1]['mean_target']/slows[0]['mean_target']-1
            check(name+' slow matched entropy',abs(slows[1]['entropy']-slows[0]['entropy']),1e-10)
        output['models'][name]=modelout
        print(name,'gain',modelout['matched_entropy_gain'],'matched amplitude',match,flush=True)
        (ROOT/'statistics_progress.json').write_text(json.dumps(output,indent=2))
    output['all_checks_pass']=all(x['pass']for x in output['checks']);output['elapsed_seconds']=time.perf_counter()-start
    (ROOT/'statistics_results.json').write_text(json.dumps(output,indent=2));print(json.dumps({'elapsed':output['elapsed_seconds'],'pass':output['all_checks_pass'],'checks':len(output['checks'])},indent=2))
    assert output['all_checks_pass']

if __name__=='__main__':main()
