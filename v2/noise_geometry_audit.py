"""Independent exact-rational audit of Cyclecraft's current-noise geometry.

No floating-point eigensolver or finite differences are used.  A frozen
reversible generator is checked by two independent constructions:
  (1) first/second tilted-eigenvalue perturbation equations;
  (2) the backward Poisson corrector and martingale quadratic variation.
Both must agree exactly with the weighted cycle-space noise floor.

Run: python v2/noise_geometry_audit.py
"""
from pathlib import Path
import json, random, time
import sympy as sp

HERE=Path(__file__).resolve().parent
CORE=[(0,1),(1,2),(2,0)]
ARCH={
    'star':CORE+[(0,3),(0,4)],
    'chain':CORE+[(0,3),(3,4)],
    'cycle':CORE+[(0,3),(3,4),(4,0)],
}

def serial(M):
    return [[str(x) for x in M.row(i)] for i in range(M.rows)]

def graph_data(edges):
    D=sp.zeros(5,len(edges))
    for e,(i,j) in enumerate(edges):D[i,e]=-1;D[j,e]=1
    C=sp.Matrix.hstack(*D.nullspace())
    return D,C

def check_case(name,edges,q,w,retain=False):
    n,m=5,len(edges);one=sp.ones(n,1)
    D,C=graph_data(edges)
    Q=sp.diag(*q);W=sp.diag(*w);L=D*Q*D.T
    P0=one*one.T/n;Lp=(L+P0).inv()-P0
    J=-Q*D.T*W;K=D*J
    Z=sum(1/x for x in w);pi=sp.Matrix([1/(Z*x) for x in w])
    normalization=pi*one.T
    group=(K+normalization).inv()-normalization
    A=sp.diag(*[2*x/Z for x in q])
    # Differentiate a signed count on each individual oriented edge.
    first=[];second=[];response=[];U=[]
    for e,(i,j) in enumerate(edges):
        R=sp.zeros(n);S=sp.zeros(n)
        R[j,i]=q[e]*w[i];R[i,j]=-q[e]*w[j]
        S[j,i]=q[e]*w[i];S[i,j]=q[e]*w[j]
        first.append(R);second.append(S)
        # Eigenvector perturbation: K r_e = -K'_e pi, 1^T r_e=0.
        response.append(-group*R*pi)
        b=R.T*one
        # Frozen corrector K^T u=-b, pi^T u=0.
        U.append(-(K.T+one*pi.T).inv()*b)
    V=sp.zeros(m)
    for a in range(m):
        for b in range(m):
            diagonal=(one.T*second[a]*pi)[0] if a==b else 0
            V[a,b]=diagonal+(one.T*first[a]*response[b])[0]+(one.T*first[b]*response[a])[0]
    U=sp.Matrix.hstack(*U)
    B=sp.eye(m)+D.T*U
    martingale=sp.simplify(B.T*A*B)
    floor=sp.simplify(C*(C.T*A.inv()*C).inv()*C.T)
    electrical=sp.simplify(sp.Rational(2,1)/Z*(Q-Q*D.T*Lp*D*Q))
    tests={
        'column_generator':one.T*K==sp.zeros(1,n),
        'equilibrium':K*pi==sp.zeros(n,1),
        'zero_frozen_mean':all((one.T*R*pi)[0]==0 for R in first),
        'cycle_basis_closed':D*C==sp.zeros(n,C.cols),
        'tilted_poisson_agree':sp.simplify(V-martingale)==sp.zeros(m),
        'tilted_graph_floor_agree':sp.simplify(V-floor)==sp.zeros(m),
        'electrical_floor_agree':sp.simplify(V-electrical)==sp.zeros(m),
        'covariance_symmetric':V==V.T,
        'covariance_rank_equals_cycle_rank':V.rank()==m-n+1,
        'gradient_kernel':sp.simplify(V*D.T)==sp.zeros(m,n),
        'cycle_covariance_positive':all((C.T*A.inv()*C)[:k,:k].det()>0 for k in range(1,C.cols+1)),
        'poisson_equations':sp.simplify(K.T*U+sp.Matrix.hstack(*[R.T*one for R in first]))==sp.zeros(n,m),
    }
    bridge_checks=[]
    for e,(tail,head) in enumerate(edges):
        adjacency={i:[] for i in range(n)}
        for f,(i,j) in enumerate(edges):
            if f!=e:adjacency[i].append(j);adjacency[j].append(i)
        seen={head};stack=[head]
        while stack:
            for j in adjacency[stack.pop()]:
                if j not in seen:seen.add(j);stack.append(j)
        if tail not in seen:
            indicator=sp.Matrix([int(i in seen) for i in range(n)])
            expected=sp.zeros(m,1);expected[e]=1
            bridge_checks.append({'edge':e,'component':sorted(seen),
                                  'gradient_equals_edge_count':D.T*indicator==expected,
                                  'zero_variance_exact':V[e,e]==0})
    tests['bridge_pathwise_potentials']=all(a['gradient_equals_edge_count'] for a in bridge_checks)
    tests['bridge_zero_diffusion']=all(a['zero_variance_exact'] for a in bridge_checks)
    out={'case':name,'architecture':name.split('_')[0],'edges':edges,'q':[str(x)for x in q],
         'w':[str(x)for x in w],'checks':tests,'bridge_checks':bridge_checks,
         'rank':V.rank(),'cycle_rank':m-n+1,'pass':all(tests.values())}
    if retain:
        out.update({'incidence':serial(D),'cycle_basis':serial(C),'equilibrium':serial(pi),
                    'edge_covariance_per_unit_time':serial(V),'noise_floor':serial(floor),
                    'poisson_correctors':serial(U),'corrected_jump_cochains':serial(B)})
    return out

def main():
    start=time.perf_counter();rng=random.Random(20261031);cases=[]
    for architecture,edges in ARCH.items():
        m=len(edges)
        cases.append(check_case(architecture+'_unit',edges,[sp.S.One]*m,[sp.S.One]*5,True))
        for k in range(4):
            q=[sp.Rational(rng.randrange(1,8),rng.randrange(1,6))for _ in edges]
            w=[sp.Rational(rng.randrange(1,7),rng.randrange(1,5))for _ in range(5)]
            cases.append(check_case(architecture+'_rational_'+str(k),edges,q,w,k==0))
    proof={
        'scope':'Connected finite fixed graph; positive periodic rates and periodic probability state; time-independent antisymmetric jump weights.',
        'periodic_corrector':'(partial_t+K^T)u_c=(mu_c/T)1-b_c, u_c(T)=u_c(0), b_c(i)=sum_j k_{j<-i} c_{i->j}.',
        'martingale':'N_c+u_c(X_t,t)-u_c(X_0,0)-(mu_c/T)t is a martingale.',
        'variance_per_period':'V_c=int_0^T sum_e a_e[c_e+(D^T u_c)_e]^2 dt, a_e=p_tail k_forward+p_head k_reverse.',
        'matrix_floor':'V_edge >= int C(C^T diag(1/a) C)^(-1)C^T dt in Loewner order.',
        'floor_remainder':'For phi*=-(D A D^T)^+D A c, V_c-floor_c=int(u_c-phi*)^T D A D^T(u_c-phi*)dt >=0.',
        'kernel':'ker V_edge=im D^T, rank V_edge=m-n+1; follows from positivity of every a_e and the bounded endpoint identity for gradients.',
        'frozen_equality':'At frozen detailed balance a=2q/Z, the Poisson corrector attains the instantaneous weighted least-squares minimum; all three exact matrices coincide.',
        'bridge_identity':'For bridge e and head component S_e, N_e(t)=1_{X_t in S_e}-1_{X_0 in S_e}. Thus |N_e|<=1 for all t; its asymptotic variance and covariance rates vanish.',
        'aggregate_identity':'If a storage subgraph is attached through one core vertex, aggregating all storage mass at that vertex gives j_core=H_core dot(p_aggregate); full statistics still use the full graph.',
        'shortcut_hazard':'If total outgoing shortcut rate <=epsilon in every state and time, killed survival S obeys S_prime>=-epsilon S. Therefore P(any shortcut by H)<=1-exp(-epsilon H)<=epsilon H. Every new auxiliary winding must use a shortcut.',
        'novelty_boundary':'Stationary current gauge invariance is established, e.g. Wachtel, Vollmer, Altaner, PRE92 042132 (2015). Here the explicit periodic floor and its architecture use are derived, not asserted globally novel.'
    }
    out={'schema':'cyclecraft-v2-noise-geometry-v1','status':'PASS'if all(x['pass']for x in cases)else'FAIL',
         'exact_cases':len(cases),'exact_checks':sum(len(x['checks'])for x in cases),'cases':cases,'proof_notes':proof,
         'elapsed_seconds':time.perf_counter()-start}
    (HERE/'noise_geometry_results.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({'status':out['status'],'cases':out['exact_cases'],'checks':out['exact_checks'],
                      'unit_ranks':{c['architecture']:c['rank']for c in cases if c['case'].endswith('_unit')},
                      'seconds':out['elapsed_seconds']},indent=2))
    assert out['status']=='PASS'

if __name__=='__main__':main()
