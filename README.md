# TROMSØ

**Tromsø Reference Orbital Model for Space Objects**

A Rust model of small-debris production, random fragment ejection, atmospheric
survival and orbital transport, conditioned on a cached Space-Track catalogue.
Built in Tromsø, Norway, for space-object radar and space-hazard studies.

## Full LEO domain

- Reporting altitude: **0--2000km**, normally25km bins, with an exact120km
  boundary. Every product includes the full coordinate array.
- Orbital survival is modeled above120km. Lower cells are NaN with a validity
  mask: they are a reentry/ablation domain, not measured zero density. Ablation
  and ground-impact modeling are not included in this initial release.
- Source parents: bound catalogue orbits intersecting LEO, including eccentric
  orbits with apogees above2000km. Current GP geometry is a stationary source
  prior, not a reconstructed collision history.
- Bound fragment trajectories are **not discarded at2000or2500km apogee**.
  They may leave and reenter LEO. DOP853 includes central+J2 Earth gravity,
  Sun/Moon third bodies from SOFA epv00/moon98, and corotating-atmosphere drag.
- Drag is neglected above2500km; this does not remove the trajectory. Solar
  radiation pressure is not yet included, so high-A/m long-term evolution
  remains a sensitivity calculation.
- Six logarithmic size bins from1mmto79mm, with the empirical small-fragment
  area-to-mass model as default. Compact spheres are a comparison, not a
  blanket morphology assumption. Characteristic length is not an RCS diameter.
- Production normalization is adjustable. Default one fragment/day is a unit
  normalization, **not an estimate of the real global debris source rate**.

## Quick start

```sh
git clone git@github.com:jvierine/tromso.git
cd tromso
cargo test --lib
cargo run --release --bin tromso -- \
  --catalogue /path/to/cached-gp.h5 --output output/full-leo-30d.h5 \
  --particles 384 --history-days 30 --threads 4
cargo run --release --bin plot_leo -- output/full-leo-30d.h5
```

No Space-Track credentials or provider traffic are used. Supply an existing
HDF5 GP snapshot containing `records_json`; snapshots and generated products
are excluded from Git. All source and scientific processing is Rust; figures
use Rustplotlib, orbit integration uses numerical_propagator, and the atmosphere
uses the native Rust rmsis implementation of NRLMSISE-00. Dependencies are pinned.

`leo/weighted_objects`, `density_per_km3_per_ln_length`, `weighted_sum_squares`
and `monte_carlo_counts` are altitude-by-size arrays. `orbital_population_valid`
is the altitude validity mask. The map leaves unsampled cells blank as well:
finite Monte Carlo zeros do not establish zero physical density. Increase the
sample count before interpreting small-scale features or estimating flux.

## Detailed physical prescription

A reproducible Rust source-and-survival model, using the existing
`numerical_propagator` DOP853 library and native `rmsis` NRLMSISE-00.
This is a sensitivity model, not a calibrated replacement for ORDEM or MASTER.
No Python, provider downloads, catalogue requests, or persistent voltages.

## Physical prescription

1. Read a cached Space-Track GP snapshot (`records_json` in HDF5). Select LEO
   parents with mean perigee120--2000km without an apogee ceiling. The default weights
   every selected object equally, including active payloads. `--parents debris`
   changes this prior. Catalogue completeness and collision-rate weighting are
   not inferred from a single snapshot.
2. Use SGP4 at each GP epoch to obtain a state, rotate TEME to GCRS, then draw
   a uniform mean anomaly on that osculating Kepler orbit. The current snapshot
   defines a stationary source geometry, **not the historical orbit of a parent
   at an invented collision date**. No SGP4 propagation backwards by years.
3. Draw characteristic lengths from a truncated power law over 1--79 mm.
   Collision cumulative exponent 1.71 (differential exponent 2.71), explosion
   cumulative exponent 1.6. Stratification allocates equal computational
   samples to six logarithmic size bins; weights correct for the actual law.
4. Draw an area-to-mass ratio, and a correlated ejection speed. The default
   `sbm-small` uses the published <8cm lognormal A/m prescription; its median
   is nearly size independent below approximately 17.8 mm. `sphere` uses
   compact aluminium spheres with density2700kg/m3 and A/m=3/(2rho_material*d).
   These are deliberately different assumptions, not interchangeable sizes.
   The NASA characteristic length is NOT a radar-equivalent diameter/RCS.
5. Ejection direction is isotropic. Collision log10(delta-v/[m/s]) has mean
   0.9 log10(A/m/[m2/kg])+2.9 and sigma0.4; explosion mean is0.2chi+1.85.
   No arbitrary cap on delta-v; unbound ejections are explicit output cases.
6. Draw source ages uniformly over the specified finite history. Propagate
   each fragment forward from birth to the evaluation date with DOP853,
   central+J2 gravity, Cd*A/m and a corotating atmosphere. Stop at120km
   geodetic altitude. Inclination, eccentricity and perigee evolve dynamically;
   no altitude-dependent small-object cutoff is imposed.
7. For constant birth rate Q, a sample in size bin b has weight
   `Q * history_days * probability(bin b) / samples(bin b)`. The sum of weights
   for surviving fragments estimates the finite-history present population
   for that Q. The default Q=1/day is a UNIT NORMALIZATION, not an estimate of
   global debris production. unbound
   ejections, numerical failures and reentries remain explicitly distinguishable.

The population transport is therefore Q followed by orbital survival and
transport, rather than an unqualified extrapolation of catalogue counts.
A snapshot already reflects large-object survival: using it as a birth prior
is an assumption and can inherit survivor bias. A physical collision source
would additionally use parent masses, pair encounter rates and event history.
No self-consistent future collisional cascade or event mass conservation is
claimed by this statistical source construction.

## Atmosphere and propagation fidelity

The initial atmosphere is a reproducible NRLMSISE-00 gtd7d density table:
5 latitude nodes uniformly spaced in sin(latitude), four UTs, two equinox
seasons, F107=F107a=150sfu and Ap=4 by default. Arithmetic mean densities are
log-interpolated every5km over120--2500km. This is a climatological sensitivity
scenario, **not observed 2026 weather or a lifetime forecast**. High-altitude
MSIS extrapolation and missing solar radiation pressure limit physical scope.
Set `--f107 70` or `--f107 220` for solar-activity sensitivity runs.

All orbit/state arithmetic is f64. Frames use IAU2006/2000A, TT-UTC69.184s,
UT1=UTC and zero polar motion, with six-hour frame tables. These approximations
are adequate for a population prior but not precise orbit determination.
no SRP, neutral winds or fragment collisions are currently included.
Returned integration errors do not establish environmental accuracy.

## Run

From this directory, with an existing private GP snapshot:

```sh
cargo test --lib
cargo run --release --bin tromso -- \
  --catalogue data/catalogue-2026-11-03.h5 \
  --output output/catalogue-sbm-365d.h5 \
  --particles 1536 --history-days 365 --threads 4
cargo run --release --bin tromso -- \
  --catalogue data/catalogue-2026-11-03.h5 \
  --output output/lifetime-sphere-365d.h5 \
  --lifetime-grid --area-mass sphere --history-days 365 --threads 4
cargo run --release --bin tromso -- \
  --catalogue data/catalogue-2026-11-03.h5 \
  --output output/lifetime-sbm-365d.h5 \
  --lifetime-grid --area-mass sbm-small --history-days 365 --threads 4
cargo run --release --bin plot
```

The lifetime grid uses geocentric initial heights150,200,300,400,600,800,1000,1500,2000km and
lengths1mm,1cm,5cm; inclination53deg and zero ejection. For `sbm-small`, it
uses the MEDIAN A/m rather than random draws, so equal-median cases are paired.
It reports exact reentry events or censored lifetimes >365days, not invented
reentry times for survivors. Triangles in lifetime plots denote lower limits.
The catalogue Monte Carlo samples the full lognormal A/m and ejection laws.

HDF5 output retains per-fragment source ID, size, weight, A/m, mass, impulse,
initial/final Cartesian state, status, event time and integrator evaluation
count, plus catalogue SHA256 and model assumptions. A/m grid and force constants
are retained. No CSV products. Source snapshots stay ignored and private.

## Checks

Six tests cover power-law decade ratios, published small-fragment A/m
boundaries, isotropic directions, orbital-element closure, and numerical drag
energy loss with120s vs60s maximum-step convergence. Population and lifetime
runs also check numerical return status and endpoint/event coverage. These are
physical/numerical checks, not external validation of a debris environment.

## Literature

- Johnson et al.(2001), NASA's new breakup model of EVOLVE4.0,
  https://doi.org/10.1016/S0273-1177(01)00423-9 .
- Apetrii et al.(2024), explicit size, A/m and ejection laws and Table7,
  https://doi.org/10.1007/s10569-024-10205-3 .
- Picone et al.(2002), NRLMSISE-00,
  https://doi.org/10.1029/2002JA009430 .
- NASA RP1408, orbital-debris lifetime, A/m and perigee dependence,
  https://ntrs.nasa.gov/api/citations/19970034583/downloads/19970034583.pdf .

The cumulative -1.71 production law is not automatically the present-day size
law at any altitude. For compact equal-density spheres, A/m scales as1/d and
circular drag residence time approximately asd; a drag-limited steady source
can therefore flatten the logarithmic-bin production slope by about one.
For the small-fragment NASA A/m prescription, the same scaling does not follow.
Downward transport from higher orbits, source history, density and perigee must
be included before asserting an altitude-dependent small/large ratio.

## Free-molecular coefficient convention

The propagator uses `a_D = -0.5 * Cd * rho * (A/m) * |v_rel| * v_rel`.
This is the aerodynamic definition of Cd, not a continuum-flow assumption.
In the cold, hypersonic, fully absorbing limit, incident molecular momentum
flux gives `F = rho * v^2 * A`, so Cd=2 under this convention. The nominal
Cd=2.2 corresponds to momentum-flux coefficient Gamma=Cd/2=1.1. Do not remove
the half without also converting Cd. Gas/surface accommodation, temperature,
shape and tumbling can change the coefficient; this model does not calculate
a shape-resolved molecular-scattering integral.

The empirical breakup A/m prescription remains the main population model;
compact spheres are a sensitivity comparison. There is no universal 1mm
shape transition. DebriSat measurements show material-dependent plates,
rods/needles and nuggets. See NASA's experiment/fragment characterization:
https://ntrs.nasa.gov/api/citations/20220009373/downloads/Cowardin_DebriSat%20Hypervelocity_Fragments_HVIS%202022_FINAL.pdf

## Initial full-domain check

A cached2026-11-03 catalogue supplied29,771eligible LEO/intersecting parents.
The initial384weighted fragments over30days completed with0numerical failures:
229retained bound survivors,150reentries and5unbound ejections. Bound survivors
outside2000km remain propagated; only their occupancy inside LEO enters the
reported grid. This small run checks the data path and domain handling, not
population convergence, calibrated flux, or agreement with ORDEM/MASTER.

Use `--epoch 2026-11-03T12:00:00Z` to set the model evaluation epoch explicitly.
The climatological atmosphere is independent of that date; exact historical
solar/geomagnetic drivers are not yet an option in this release.
