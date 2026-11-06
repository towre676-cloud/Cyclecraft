"""Exact algebra for Cyclecraft II: acyclic storage and productive circulation.

Run with Python and SymPy. Every reported identity is rational or symbolic;
there is no finite difference, numerical rank tolerance, or sampled proof.
Column incidence has -1 at the tail and +1 at the head. At the base E=0,
G=n(d pi)^T L^+ (d pi), and B_01 is tied to control E_2.
"""
from pathlib import Path
import json
import random
import time
import sympy as s


def encode(matrix):
    return [[str(v) for v in row] for row in matrix.tolist()]


def zero(matrix):
    return all(v == 0 for v in matrix)


def graph(n, edges, conductances):
    D = s.zeros(n, len(edges))
    for e, (i, j) in enumerate(edges):
        D[i, e] = -1
        D[j, e] = 1
    Q = s.diag(*conductances)
    L = D * Q * D.T
    constant = s.ones(n, n) / n
    R = (L + constant).inv() - constant
    H = Q * D.T * R
    return D, Q, L, R, H


def zero_mass_basis(n):
    T = s.zeros(n, n - 1)
    for k in range(n - 1):
        T[0, k] = -1
        T[k + 1, k] = 1
    return T


def component_without_edge(n, edges, removed, start):
    adjacency = [[] for _ in range(n)]
    for e, (i, j) in enumerate(edges):
        if e == removed:
            continue
        adjacency[i].append(j)
        adjacency[j].append(i)
    seen = {start}
    stack = [start]
    while stack:
        i = stack.pop()
        for j in adjacency[i]:
            if j not in seen:
                seen.add(j)
                stack.append(j)
    return sorted(seen)


def tree_storage_audit(name, core_n, core_edges, edges, q, attachment):
    n = 1 + max(max(e) for e in edges)
    core_m = len(core_edges)
    D, Q, L, R, H = graph(n, edges, q)
    Dc, Qc, Lc, Rc, Hc = graph(core_n, core_edges, q[:core_m])
    T = zero_mass_basis(n)
    aggregate = s.zeros(core_n, n)
    for i in range(core_n):
        aggregate[i, i] = 1
    for i in range(core_n, n):
        aggregate[attachment, i] = 1
    energy_tree = s.zeros(n, n)
    bridges = []
    for e in range(core_m, len(edges)):
        head_component = component_without_edge(n, edges, e, edges[e][1])
        f = s.Matrix([int(i in head_component) for i in range(n)])
        unit = s.eye(len(edges))[:, e]
        identity = zero(D.T * f - unit)
        energy_tree += f * f.T / q[e]
        bridges.append({"edge_index": e, "edge": list(edges[e]),
                        "head_component": head_component,
                        "potential": [int(x) for x in f],
                        "gradient_is_edge_indicator": identity})
    checks = {
        "connected_incidence_rank": D.rank() == n - 1,
        "aggregate_routing_all_zero_mass_velocities": zero(
            (H[:core_m, :] - Hc * aggregate) * T),
        "resistance_energy_split_all_zero_mass_velocities": zero(
            T.T * (R - aggregate.T * Rc * aggregate - energy_tree) * T),
        "every_auxiliary_signed_count_is_endpoint_difference": all(
            b["gradient_is_edge_indicator"] for b in bridges),
        "cycle_rank_unchanged_by_storage": len(edges) - n + 1 == core_m - core_n + 1,
        "aggregate_normalization": zero(s.ones(1, core_n) * aggregate - s.ones(1, n)),
    }
    return {"name": name, "vertices": n, "core_vertices": core_n,
            "attachment": attachment, "edges": [list(e) for e in edges],
            "conductances": [str(v) for v in q], "aggregate_matrix": encode(aggregate),
            "bridges": bridges, "checks": checks, "pass": all(checks.values())}


def local_geometry(name, edges, n, expected_base_squared, expected_full_squared):
    d = n - 1
    D, Q, L, R, H = graph(n, edges, [s.Integer(1)] * len(edges))
    U = s.zeros(n, d)
    for a in range(d):
        U[a + 1, a] = 1
    V = s.zeros(len(edges), d)
    V[0, 1] = 1
    pi_gradient = -(s.eye(n) - s.ones(n, n) / n) * U / n
    G = n * pi_gradient.T * R * pi_gradient
    cycle_projection = s.eye(len(edges)) - H * D
    H_derivatives = [-cycle_projection * s.diag(*list(V[:, a])) * H
                     for a in range(d)]
    F = []
    for e in range(len(edges)):
        Fe = s.zeros(d, d)
        for a in range(d):
            for b in range(d):
                Fe[a, b] = (H_derivatives[a] * pi_gradient[:, b]
                            - H_derivatives[b] * pi_gradient[:, a])[e]
        F.append(Fe)
    z = s.symbols("z")
    char_base = s.factor((G[:2, :2].inv() * F[0][:2, :2]).charpoly(z).as_expr())
    char_full = s.factor((G.inv() * F[0]).charpoly(z).as_expr())
    checks = {
        "metric_symmetric": G == G.T,
        "metric_positive_definite_exact_minors": all(G[:i, :i].det() > 0 for i in range(1, d + 1)),
        "all_edge_curvatures_antisymmetric": all(zero(f + f.T) for f in F),
        "active_curvature_conserved": F[0] == F[1] == F[2],
        "auxiliary_edge_curvatures_zero": all(zero(f) for f in F[3:]),
        "active_characteristic_exact": s.expand(char_base - (z*z + expected_base_squared)) == 0,
        "assisted_characteristic_exact": s.expand(
            char_full - z**(d - 2) * (z*z + expected_full_squared)) == 0,
        "strict_assistance_exact": expected_full_squared > expected_base_squared,
        "metric_identity": zero(H.T * Q.inv() * H - R),
        "finite_speed_routing_identity": zero(H * D * Q * D.T - Q * D.T),
    }
    degrees = [sum(int(i in edge) for edge in edges) for i in range(n)]
    ratio = s.sqrt(expected_full_squared / expected_base_squared)
    return {"name": name, "vertices": n, "edges": [list(e) for e in edges],
            "conductances": ["1"] * len(edges), "well_control_matrix": encode(U),
            "barrier_control_matrix": encode(V), "incidence": encode(D),
            "laplacian_pseudoinverse": encode(R), "transport_map": encode(H),
            "equilibrium_probability_gradient": encode(pi_gradient),
            "metric": encode(G), "target_curvature": encode(F[0]),
            "metric_leading_principal_minors": [str(G[:i, :i].det()) for i in range(1, d + 1)],
            "active_characteristic": str(char_base), "assisted_characteristic": str(char_full),
            "kappa_active_squared": str(expected_base_squared),
            "kappa_assisted_squared": str(expected_full_squared),
            "kappa_active": str(s.N(s.sqrt(expected_base_squared), 40)),
            "kappa_assisted": str(s.N(s.sqrt(expected_full_squared), 40)),
            "assistance_ratio_exact": str(s.simplify(ratio)),
            "assistance_gain_fraction": str(s.N(ratio - 1, 40)),
            "baseline_resource_accounting": {
                "state_count": n, "edge_count": len(edges),
                "total_installed_conductance": len(edges),
                "uniform_probability": str(s.Rational(1, n)),
                "expected_jump_activity_rate": str(s.Rational(2 * len(edges), n)),
                "state_escape_rates": degrees,
                "maximum_escape_rate": max(degrees)},
            "checks": checks, "pass": all(checks.values())}


def cyclic_auxiliary_audit():
    core_edges = [(0, 1), (1, 2), (2, 0)]
    edges = core_edges + [(0, 3), (3, 4), (4, 0)]
    q = [s.Rational(1, 3), s.Rational(2, 5), s.Rational(7, 4),
         s.Rational(3, 2), s.Rational(5, 7), s.Rational(11, 6)]
    D, Q, L, R, H = graph(5, edges, q)
    Dc, Qc, Lc, Rc, Hc = graph(3, core_edges, q[:3])
    aux_local_edges = [(0, 1), (1, 2), (2, 0)]
    Da, Qa, La, Ra, Ha = graph(3, aux_local_edges, q[3:])
    aggregate = s.Matrix([[1, 0, 0, 1, 1], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0]])
    # Local auxiliary source is (-b3-b4,b3,b4); the cut-root balances its mass.
    aux_source = s.Matrix([[0, 0, 0, -1, -1], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
    T = zero_mass_basis(5)
    checks = {
        "aggregate_routing_independent_of_internal_auxiliary_cycle": zero((H[:3, :] - Hc * aggregate) * T),
        "articulation_energy_split": zero(T.T * (R - aggregate.T * Rc * aggregate
                                                  - aux_source.T * Ra * aux_source) * T),
        "auxiliary_cycle_adds_homology_dimension": len(edges) - 5 + 1 == 2,
        "auxiliary_cycle_observable_not_gradient": (D.T).row_join(s.eye(6)[:, 3]).rank() > D.T.rank(),
    }
    return {"name": "cyclic_auxiliary_articulation", "edges": [list(e) for e in edges],
            "conductances": [str(v) for v in q], "aggregate_matrix": encode(aggregate),
            "auxiliary_balanced_source_map": encode(aux_source),
            "checks": checks, "pass": all(checks.values())}


def bridge_mgf_audit():
    """Finite-horizon identity on exact rational uniformization kernels.

    A uniformization kernel is an exact discrete transition kernel with the
    graph's edges, rather than a floating approximation to a period ODE solve.
    The same diagonal similarity proves the theorem for the actual continuous
    period propagator. The rational examples check endpoint correlations and
    all signed-count probabilities without numerical differentiation.
    """
    z, lam = s.symbols("z lam", nonzero=True)
    cases = []
    examples = [
        ("one_leaf", 4, [(0, 1), (1, 2), (2, 0), (0, 3)]),
        ("two_leaves", 5, [(0, 1), (1, 2), (2, 0), (0, 3), (0, 4)]),
        ("chain", 5, [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4)]),
    ]
    for name, n, edges in examples:
        D, Q, L, R, H = graph(n, edges, [s.Integer(1)] * len(edges))
        max_degree = max(L.diagonal())
        M = s.eye(n) - L / (2 * max_degree)
        pi = s.ones(n, 1) / n
        transient = s.Matrix([s.Rational(i + 1, n * (n + 1) // 2) for i in range(n)])
        for e in range(3, len(edges)):
            component = component_without_edge(n, edges, e, edges[e][1])
            f = [int(i in component) for i in range(n)]
            Z = s.diag(*[z ** f_i for f_i in f])
            tilted = s.Matrix(n, n, lambda j, i: M[j, i] * z ** (f[j] - f[i]))
            checks = {
                "kernel_stochastic_exact": zero(s.ones(1, n) * M - s.ones(1, n)),
                "stationary_probability_exact": zero(M * pi - pi),
                "tilted_kernel_diagonal_similarity": zero(tilted - Z * M * Z.inv()),
                "tilted_characteristic_independent_of_count_field": s.expand(
                    tilted.charpoly(lam).as_expr() - M.charpoly(lam).as_expr()) == 0,
                "gradient_potential_is_single_bridge_count":
                    zero(D.T * s.Matrix(f) - s.eye(len(edges))[:, e]),
            }
            endpoint_examples = []
            for steps in [1, 3, 12]:
                Mn = M ** steps
                for label, p0 in [("periodic", pi), ("transient", transient)]:
                    plus = sum(Mn[j, i] * p0[i] for i in range(n) for j in range(n)
                               if f[i] == 0 and f[j] == 1)
                    minus = sum(Mn[j, i] * p0[i] for i in range(n) for j in range(n)
                                if f[i] == 1 and f[j] == 0)
                    stay = 1 - plus - minus
                    mgf = (s.ones(1, n) * Z * Mn * Z.inv() * p0)[0]
                    polynomial = stay + plus * z + minus / z
                    mean = plus - minus
                    variance = plus + minus - mean ** 2
                    # Euler operator z d/dz is the derivative in chi=log z.
                    derivative = polynomial
                    moment_checks = []
                    for order in range(1, 7):
                        derivative = z * s.diff(derivative, z)
                        expected = plus + ((-1) ** order) * minus
                        moment_checks.append(s.simplify(derivative.subs(z, 1) - expected) == 0)
                    key = f"{label}_{steps}"
                    checks[key + "_mgf_probability_identity"] = s.simplify(mgf - polynomial) == 0
                    checks[key + "_nonnegative_probability_mass"] = min(plus, minus, stay) >= 0
                    checks[key + "_six_raw_moments_exact"] = all(moment_checks)
                    checks[key + "_variance_bounded_by_one"] = 0 <= variance <= 1
                    if label == "periodic":
                        checks[key + "_zero_mean_with_correlated_endpoints"] = plus == minus
                        checks[key + "_variance_formula"] = variance == 2 * plus
                    endpoint_examples.append({"steps": steps, "initial_condition": label,
                                              "probability_minus_one": str(minus),
                                              "probability_zero": str(stay),
                                              "probability_plus_one": str(plus),
                                              "mgf_laurent_polynomial": str(polynomial),
                                              "mean": str(mean), "variance": str(variance)})
            mass = sum(pi[i] for i in component)
            cases.append({"name": name, "edge_index": e, "component": component,
                          "uniformization_kernel": encode(M), "stationary_mass": str(mass),
                          "mixing_limit_variance": str(2 * mass * (1 - mass)),
                          "endpoint_examples": endpoint_examples,
                          "checks": checks, "pass": all(checks.values())})
    return {"theorem": "G_n(z)=1^T Z(z) M^n Z(z)^-1 p0=r0+r_plus*z+r_minus/z",
            "periodic_state": "r_plus=r_minus; finite variance=2*r_plus<=1; all odd raw moments vanish",
            "asymptotic_cumulants": "M_z is diagonally similar to M; Perron root is identically one; all asymptotic cumulants vanish",
            "mixing_limit": "variance tends to 2*p(S)*(1-p(S)), not to zero",
            "rational_audit_scope": "Exact uniformization transition kernels on each tree graph; theorem also holds for continuous-time period propagators.",
            "cases": cases, "pass": all(c["pass"] for c in cases)}


def main():
    started = time.perf_counter()
    special = [
        local_geometry("one_leaf_four_states", [(0, 1), (1, 2), (2, 0), (0, 3)], 4,
                       s.Rational(25, 162), s.Rational(11, 27)),
        local_geometry("two_leaves_five_states", [(0, 1), (1, 2), (2, 0), (0, 3), (0, 4)], 5,
                       s.Rational(49, 351), s.Rational(13, 27)),
        local_geometry("chain_five_states", [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4)], 5,
                       s.Rational(49, 513), s.Rational(11, 27)),
    ]
    structural = []
    for name, aux in [("star_split", [(0, 3), (0, 4)]),
                      ("chain_split", [(0, 3), (3, 4)])]:
        core = [(0, 1), (1, 2), (2, 0)]
        structural.append(tree_storage_audit(name, 3, core, core + aux,
                          [s.Rational(1, 5), s.Rational(2, 3), s.Rational(7, 4),
                           s.Rational(3, 2), s.Rational(11, 7)], 0))
    rng = random.Random(20261030)
    for case in range(20):
        nc = 3 + case % 3
        na = 1 + case % 4
        core = [(i, (i + 1) % nc) for i in range(nc)]
        if nc >= 4 and case % 2:
            core.append((0, 2))
        attach = case % nc
        edges = list(core)
        for j in range(na):
            node = nc + j
            parent = attach if j == 0 else rng.choice([attach] + list(range(nc, node)))
            edges.append((parent, node))
        q = [s.Rational(rng.randrange(1, 8), rng.randrange(1, 6)) for _ in edges]
        structural.append(tree_storage_audit(f"rational_case_{case:02d}", nc, core, edges, q, attach))
    cyclic = cyclic_auxiliary_audit()
    bridge_mgf = bridge_mgf_audit()
    checks = {"three_exact_local_models_pass": all(x["pass"] for x in special),
              "twenty_random_rational_graphs_pass": all(x["pass"] for x in structural[2:]),
              "star_chain_split_pass": all(x["pass"] for x in structural[:2]),
              "cyclic_auxiliary_split_pass": cyclic["pass"],
              "finite_horizon_bridge_mgf_pass": bridge_mgf["pass"],
              "star_chain_same_states_edges_installed_conductance_activity":
              all(special[1]["baseline_resource_accounting"][k] == special[2]["baseline_resource_accounting"][k]
                  for k in ["state_count", "edge_count", "total_installed_conductance", "expected_jump_activity_rate"])}
    checks_count = (sum(len(x["checks"]) for x in special + structural + [cyclic])
                    + sum(len(c["checks"]) for c in bridge_mgf["cases"]) + len(checks))
    out = {"schema": "cyclecraft.v2.exact-tree-audit.v1", "status": "PASS" if all(checks.values()) else "FAIL",
           "precision": "Exact SymPy rational arithmetic; decimal values are display-only 40-digit evaluations.",
           "control_convention": "E=(0,x,y,z[,w]); B_01=y; all other barriers fixed; q_e=1 and E=0 at the local base.",
           "theorems": {
               "aggregate_routing": "j_core=H_core(q_core) d(p_core with all single-port storage occupation assigned to its attachment)/dt",
               "tree_resistance_split": "b^T L^+ b=bar(b)^T L_core^+ bar(b)+sum_tree_edges b(S_e)^2/q_e, for 1^T b=0",
               "bridge_count": "N_e(t)=1{X_t in S_e}-1{X_0 in S_e}, for edge oriented into S_e; each count belongs to {-1,0,1}",
               "covariance_consequence": "Every auxiliary bridge count has bounded variance and zero long-time variance and covariance rate.",
               "scope": "Local curvature/action gains do not assert finite-amplitude or rate/activity-matched global optimality."},
           "local_models": special, "structural_audits": structural, "cyclic_auxiliary_audit": cyclic,
           "bridge_finite_horizon_mgf": bridge_mgf,
           "checks": checks, "exact_predicate_count": checks_count,
           "elapsed_seconds": time.perf_counter() - started}
    Path(__file__).with_name("tree_exact_results.json").write_text(json.dumps(out, indent=2, default=lambda x: bool(x) if isinstance(x, s.logic.boolalg.BooleanAtom) else str(x)) + "\n")
    print(json.dumps({"status": out["status"], "exact_predicate_count": checks_count,
                      "checks": checks, "elapsed_seconds": out["elapsed_seconds"],
                      "local_gains": [{"name": m["name"], "ratio": m["assistance_ratio_exact"],
                                       "gain": m["assistance_gain_fraction"]} for m in special]}, indent=2))
    if out["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
