# Cyclecraft II — Storage without winding

The companion `Cyclecraft.html` is an offline scientific paper with embedded
equations, fonts, data, autonomous figures and a downloadable numerical program.
`Cyclecraft_Complete.zip` contains the full source and the first edition archive.

The finite-speed two-leaf construction has exactly certified mean current
0.07369813425630575 per cycle. A continuous-time upper bound for every
frozen-storage comparator in the specified fixed-gate class is
0.06423038232315985: at least 14.74% advantage. Chain and cyclic auxiliary
constructions have their own certified bounds. Tree-edge signed counts are
bounded pathwise, independently of observation time.

## Replay without optimization

```
python -m pip install -r requirements-core.txt
python replay_v2.py
```

The retained optimization proposals are untrusted inputs. The replay recomputes
the exact tree/noise audits, interval probability enclosures, exact continuous
assisted feasibility, continuous baseline dual bounds, frontier bounds and
independent counting dynamics. It writes `replay_v2_receipt.json`. Default
verification takes about two minutes in the reference environment.

To regenerate the proposals and the harmonic comparison program as well:

```
python -m pip install -r requirements-opt.txt
python replay_v2.py --full
```

The proposal solver is CVXPY/Clarabel. The verifier uses exact rational algebra,
Arb enclosures, linear-program dual residual bounds and independent dynamics.
Optimizer success is not the certificate.

## Experiment families

1. Exact local geometry: unit conductances, rational metric and curvature.
2. Smooth harmonic comparisons: period 20/80, entropy-matched architectures;
   their activity budgets are measured, not silently matched.
3. Continuous resource certificates: period 8, entropy≤3/20, activity≤16,
   energies within±6/5, escape≤6, probability≥3/200, mean inventories1/5,
   assisted well slew≤20; fixed positive piecewise-linear gates.
4. Explicit dyadic three-hold protocols: period 6, rigorous finite-horizon
   counting risk and 320-bit Arb enclosures, independent 110-digit reference.

The fourth family illustrates counting protection; it is not an entropy- or
activity-matched comparison with the optimized continuous waveforms.

## Numerical sources

- `tree_exact.py`: exact geometry, aggregation, tree energy, bridge potentials,
  endpoint-correlated finite-count generating functions.
- `noise_geometry_audit.py`: exact frozen-equilibrium covariance and graph floor.
- `statistics.py`: periodic dynamics, full joint counting derivatives, backward
  Poisson corrector, finite-horizon moments, killed count lift and seeded paths.
- `precision_cert.py`: exact rational held generators, Arb exponentials, periodic
  initialization error, finite-horizon cumulants/tails and rare-edge risks.
- `convex_compile.py`: optimization proposals; CVXPY is needed only to solve anew.
- `primal_certificate.py`: exact rational repair, polynomial dynamics,
  Bernstein feasibility, rigorous entropy and slew bounds.
- `dual_certificate.py`: continuous-time upper certificate using interval
  coefficients, entropy supporting planes and certified LP-dual corrections.
- `aggregate_convex.py`: source-bound gains and resource-frontier certificates.
- `noise_smooth_review.py`: independent electrical-floor cross-check.

All generated JSON records retain source conventions, budgets and precision.
`SHA256SUMS.json` binds the packaged sources and retained evidence.

## Reading the guarantees

The global comparison is conditional on its fixed gate waveform, mean inventory
and declared resource class. The baseline permits every active well control and
freezes auxiliary wells. Its upper bound permits unlimited well slew, making it
a stronger comparator. The assisted construction satisfies finite well slew.

Tree count bounds concern signed edge crossings, not absence of visits or gross
traffic. Productive-current noise remains. All finite-time statistics use the
full network; aggregate mean routing is not an assumption of Markov lumpability.

The continuous witness and dual bounds use exact rational/interval certificates.
Independent smooth counting integrations are binary64 checks. The dyadic held
protocol probabilities have rigorous ball enclosures. These scopes are stated
separately in the paper and JSON.

Time is dimensionless and energies are in bath thermal units. Physical gate
driver expenditure and external load work require a calibrated device model.
No hardware performance measurement is asserted.

Reference environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0,
mpmath 1.3.0, python-flint 0.9.0; optional CVXPY 1.9.3, Clarabel 0.11.1.
