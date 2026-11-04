//! Catalogue-conditioned, weighted fragment-source model. Not ORDEM or MASTER.
//! Laws: Johnson et al. (2001), as tabulated in Apetrii et al. (2024),
//! DOI10.1007/s10569-024-10205-3, equations1,9-11,19-22 and Table7.
use anyhow::{Result, ensure};
use numerical_propagator::{Inputs, Stats, engine_free, engine_new, engine_propagate};
use serde::{Deserialize, Serialize};
use std::sync::OnceLock;
pub const MU: f64 = 3.986004418e14;
pub const RE: f64 = 6378137.;
pub const DAY: f64 = 86400.;
pub const YEAR: f64 = 365.25 * DAY;
pub type V = [f64; 3];
pub fn dot(a: V, b: V) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
pub fn norm(a: V) -> f64 {
    dot(a, a).sqrt()
}
pub fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub struct Rng(pub u64);
impl Rng {
    pub fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        (((z ^ (z >> 31)) >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    pub fn normal(&mut self) -> f64 {
        (-2. * self.uniform().ln()).sqrt() * (std::f64::consts::TAU * self.uniform()).cos()
    }
    pub fn direction(&mut self) -> V {
        let z = 2. * self.uniform() - 1.;
        let p = std::f64::consts::TAU * self.uniform();
        let r = (1. - z * z).sqrt();
        [r * p.cos(), r * p.sin(), z]
    }
}
pub fn power_quantile(u: f64, lo: f64, hi: f64, p: f64) -> f64 {
    ((1. - u) * lo.powf(-p) + u * hi.powf(-p)).powf(-1. / p)
}
pub fn bin_probability(lo: f64, hi: f64, min: f64, max: f64, p: f64) -> f64 {
    (lo.powf(-p) - hi.powf(-p)) / (min.powf(-p) - max.powf(-p))
}
/// Small-fragment published area/mass log10 distribution, ONLY Lc<0.08m.
pub fn small_am_parameters(d: f64) -> Result<(f64, f64)> {
    ensure!(
        d > 0. && d < 0.08,
        "Small-fragment prescription requires 0<Lc<8cm"
    );
    let l = d.log10();
    let mean = if l <= -1.75 {
        -0.3
    } else if l < -1.25 {
        -0.3 - 1.4 * (l + 1.75)
    } else {
        -1.
    };
    let sigma = if l <= -3.5 {
        0.2
    } else {
        0.2 + 0.1333 * (l + 3.5)
    };
    Ok((mean, sigma))
}
pub fn area(d: f64) -> f64 {
    if d < 0.00167 {
        0.540424 * d * d
    } else {
        0.556945 * d.powf(2.0047077)
    }
}
pub fn sphere_am(d: f64, density: f64) -> f64 {
    1.5 / (density * d)
}
pub fn sample_am(d: f64, law: &str, rng: &mut Rng) -> Result<f64> {
    match law {
        "sbm-small" => {
            let (m, s) = small_am_parameters(d)?;
            Ok(10f64.powf(m + s * rng.normal()))
        }
        "sphere" => Ok(sphere_am(d, 2700.)),
        _ => anyhow::bail!("Unknown area/mass model"),
    }
}
pub fn ejection(am: f64, collision: bool, rng: &mut Rng) -> V {
    let (a, b) = if collision { (0.9, 2.9) } else { (0.2, 1.85) };
    let speed = 10f64.powf(a * am.log10() + b + 0.4 * rng.normal());
    rng.direction().map(|x| x * speed)
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Parent {
    pub id: u64,
    pub name: String,
    pub epoch: f64,
    pub state: [f64; 6],
    pub object_type: String,
}
fn mv(a: [[f64; 3]; 3], v: V) -> V {
    a.map(|r| dot(r, v))
}
fn mt(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}
fn mm(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}
fn teme_rotation(t: f64) -> [[f64; 3]; 3] {
    let jd = t / DAY;
    let g = sofars::erst::gmst82(2440587.5, jd);
    let (c, s) = (g.cos(), g.sin());
    let r = [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]];
    let c2t = sofars::pnp::c2t06a(2440587.5, jd + 69.184 / DAY, 2440587.5, jd, 0., 0.);
    mm(mt(c2t), r)
}
/// TEME -> GCRS with UT1=UTC, zero polar motion, TT-UTC=69.184s.
/// Numerical frame-rate term included. Suitable for population prior, not POD.
pub fn teme_gcrs(t: f64, p: V, v: V) -> [f64; 6] {
    let m = teme_rotation(t);
    let p = p.map(|x| 1000. * x);
    let v = v.map(|x| 1000. * x);
    let r = mv(m, p);
    let w = mv(m, v);
    let l = mv(teme_rotation(t - 0.5), p);
    let h = mv(teme_rotation(t + 0.5), p);
    [
        r[0],
        r[1],
        r[2],
        w[0] + h[0] - l[0],
        w[1] + h[1] - l[1],
        w[2] + h[2] - l[2],
    ]
}
pub fn orbital_bounds(y: [f64; 6]) -> (f64, f64) {
    let r = [y[0], y[1], y[2]];
    let v = [y[3], y[4], y[5]];
    let rr = norm(r);
    let energy = dot(v, v) / 2. - MU / rr;
    if energy >= 0. {
        return (f64::NEG_INFINITY, f64::INFINITY);
    }
    let a = -MU / (2. * energy);
    let h = cross(r, v);
    let e = (1. + 2. * energy * dot(h, h) / (MU * MU)).max(0.).sqrt();
    (a * (1. - e) - RE, a * (1. + e) - RE)
}
#[derive(Clone)]
pub struct Atmosphere {
    pub altitude: Vec<f64>,
    pub log_density: Vec<f64>,
    pub f107: f64,
}
static ATM: OnceLock<Atmosphere> = OnceLock::new();
impl Atmosphere {
    pub fn build(f107: f64) -> Self {
        let model = rmsis::Model::new();
        let mut altitude = vec![];
        let mut log_density = vec![];
        for k in 0..=476 {
            let h = 120. + 5. * k as f64;
            let mut sum = 0.;
            let mut n = 0.;
            for day in [80, 264] {
                for sinlat in [-0.8f64, -0.4, 0., 0.4, 0.8] {
                    for ut in [0., 21600., 43200., 64800.] {
                        let input = rmsis::Input::consistent(
                            day,
                            ut,
                            h,
                            sinlat.asin().to_degrees(),
                            0.,
                            f107,
                            f107,
                            rmsis::MagneticActivity::Daily { ap: 4. },
                        );
                        let q = model.gtd7d(&input).total_mass;
                        assert!(q.is_finite() && q > 0.);
                        sum += q;
                        n += 1.;
                    }
                }
            }
            altitude.push(h);
            log_density.push((sum / n).ln());
        }
        Self {
            altitude,
            log_density,
            f107,
        }
    }
    pub fn density(&self, h: f64) -> f64 {
        // Drag above 2500km is neglected, not used to discard an orbit.
        if h > 2500. {
            return 0.;
        }
        let x = ((h - 120.) / 5.).clamp(0., (self.altitude.len() - 1) as f64);
        let i = (x.floor() as usize).min(self.altitude.len() - 2);
        let u = x - i as f64;
        ((1. - u) * self.log_density[i] + u * self.log_density[i + 1]).exp()
    }
    pub fn install(self) -> Result<()> {
        ATM.set(self)
            .map_err(|_| anyhow::anyhow!("Atmosphere already installed"))
    }
}
fn density_callback(_: f64, _: f64, _: f64, _: f64, h: f64, _: f64, _: f64, _: [f64; 7]) -> f64 {
    ATM.get().expect("Atmosphere installed").density(h)
}
#[derive(Clone)]
pub struct Environment {
    pub start: f64,
    pub end: f64,
    pub frames: Vec<f64>,
    pub polar: Vec<f64>,
    pub step: f64,
    pub f107: f64,
    pub sun_pos: Vec<f64>,
    pub sun_vel: Vec<f64>,
    pub moon_pos: Vec<f64>,
    pub moon_vel: Vec<f64>,
}
impl Environment {
    pub fn new(start: f64, end: f64, f107: f64) -> Self {
        let step = 21600.;
        let count = ((end - start) / step).ceil() as usize + 3;
        let mut frames = Vec::new();
        let mut polar = Vec::new();
        let (mut sun_pos, mut sun_vel, mut moon_pos, mut moon_vel) =
            (vec![], vec![], vec![], vec![]);
        for k in 0..count {
            let jd = (start + k as f64 * step + 69.184) / DAY;
            let q = sofars::pnp::c2i06a(2440587.5, jd);
            frames.extend(q.into_iter().flatten());
            polar.extend([1., 0., 0., 0., 1., 0., 0., 0., 1.]);
            let (earth, _) = sofars::eph::epv00(2440587.5, jd).expect("SOFA epoch range");
            let moon = sofars::eph::moon98(2440587.5, jd);
            const AU: f64 = 149597870700.;
            sun_pos.extend(earth[0].map(|x| -AU * x));
            sun_vel.extend(earth[1].map(|x| -AU * x / DAY));
            moon_pos.extend(moon[0].map(|x| AU * x));
            moon_vel.extend(moon[1].map(|x| AU * x / DAY));
        }
        Self {
            start,
            end,
            frames,
            polar,
            step,
            f107,
            sun_pos,
            sun_vel,
            moon_pos,
            moon_vel,
        }
    }
    pub fn propagate(
        &self,
        start: f64,
        end: f64,
        y: [f64; 6],
        beta: f64,
        rtol: f64,
        max_step: f64,
        sample_step: f64,
    ) -> Result<Track> {
        ensure!(
            start >= self.start && end <= self.end && end > start,
            "Time bounds"
        );
        let mut c = [0f64; 9];
        c[0] = 1.;
        c[6] = -1.08262668e-3 / 5f64.sqrt();
        let zeros = [0.; 7];
        let flux = [self.f107];
        let doy = [307.];
        let d = vec![0.; self.frames.len() / 9];
        let cs = [0.; 9];
        let input = Inputs {
            degree: 2,
            gm: MU,
            radius: RE,
            coefficients_c: c.as_ptr(),
            coefficients_s: cs.as_ptr(),
            frame_start: self.start,
            frame_step: self.step,
            nframe: self.frames.len() / 9,
            q: self.frames.as_ptr(),
            polar: self.polar.as_ptr(),
            dut: d.as_ptr(),
            body_start: self.start,
            body_step: self.step,
            nbody: self.sun_pos.len() / 3,
            sun_pos: self.sun_pos.as_ptr(),
            sun_vel: self.sun_vel.as_ptr(),
            moon_pos: self.moon_pos.as_ptr(),
            moon_vel: self.moon_vel.as_ptr(),
            weather_start: self.start,
            weather_step: self.end - self.start + DAY,
            nweather: 1,
            f107: flux.as_ptr(),
            f107a: flux.as_ptr(),
            ap: zeros.as_ptr(),
            doy: doy.as_ptr(),
            correction_start: 0.,
            correction_step: 1.,
            ncorrection: 0,
            log_correction: zeros.as_ptr(),
            correction_height_linear: zeros.as_ptr(),
            correction_height_quadratic: zeros.as_ptr(),
            correction_height_origin: 0.,
            correction_height_scale: 1.,
            correction_clip: 0.,
            beta,
            drag: (beta > 0.) as i32,
            third_bodies: 1,
        };
        let count = ((end - start) / sample_step).ceil() as usize;
        ensure!(count < 200000, "Sample cap");
        let mut times: Vec<_> = (0..count).map(|k| start + k as f64 * sample_step).collect();
        times.push(end);
        let mut t = vec![0.; times.len() + 1];
        let mut states = vec![0.; 6 * (times.len() + 1)];
        let mut stats = Stats::default();
        let code = unsafe {
            let ctx = engine_new(&input);
            ensure!(!ctx.is_null(), "Numerical context creation");
            (*ctx).set_density_callback(density_callback);
            let code = engine_propagate(
                ctx,
                start,
                end,
                y.as_ptr(),
                6,
                rtol,
                max_step,
                120.,
                times.as_ptr(),
                times.len(),
                t.as_mut_ptr(),
                states.as_mut_ptr(),
                &mut stats,
            );
            engine_free(ctx);
            code
        };
        ensure!(code == 0, "DOP853 error {code}");
        ensure!(stats.nout > 0, "No numerical output");
        t.truncate(stats.nout);
        states.truncate(6 * stats.nout);
        if !stats.event.is_finite() {
            ensure!(
                (t.last().unwrap() - end).abs() < 0.01,
                "Incomplete numerical output"
            );
        }
        Ok(Track {
            time: t,
            state: states,
            event: stats.event,
            nfev: stats.nfev,
        })
    }
}
pub struct Track {
    pub time: Vec<f64>,
    pub state: Vec<f64>,
    pub event: f64,
    pub nfev: usize,
}
/// Full LEO reporting domain. Below120km is a reentry regime, not a modeled
/// steady orbital population. Use NaN plus a validity mask, never measured zero.
pub fn leo_altitude_edges() -> Vec<f64> {
    let mut e: Vec<_> = (0..=80)
        .map(|i| i as f64 * 25.)
        .filter(|h| *h != 125.)
        .collect();
    e.push(120.);
    e.sort_by(f64::total_cmp);
    e
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn powerlaw_decade_counts_and_quantiles() {
        assert!(
            (bin_probability(0.001, 0.01, 0.001, 0.1, 1.71)
                / bin_probability(0.01, 0.1, 0.001, 0.1, 1.71)
                - 10f64.powf(1.71))
            .abs()
                < 1e-10
        );
        assert!((power_quantile(0., 0.001, 0.08, 1.71) - 0.001).abs() < 1e-14);
        assert!((power_quantile(1., 0.001, 0.08, 1.71) - 0.08).abs() < 1e-12);
    }
    #[test]
    fn published_small_am_and_sphere_differ() {
        assert_eq!(small_am_parameters(0.001).unwrap().0, -0.3);
        assert_eq!(small_am_parameters(0.01).unwrap().0, -0.3);
        assert_eq!(small_am_parameters(0.07).unwrap().0, -1.);
        assert!(small_am_parameters(0.1).is_err());
        assert!((sphere_am(0.001, 2700.) / sphere_am(0.01, 2700.) - 10.).abs() < 1e-12);
    }
    #[test]
    fn isotropic_ejection_has_no_preferred_direction() {
        let mut rng = Rng(42);
        let mut m = [0.; 3];
        for _ in 0..50000 {
            let v = rng.direction();
            assert!((norm(v) - 1.).abs() < 1e-12);
            for j in 0..3 {
                m[j] += v[j] / 50000.;
            }
        }
        assert!(m.iter().all(|x| x.abs() < 0.015));
    }
    #[test]
    fn zero_kick_preserves_bounds() {
        let r = RE + 500000.;
        let y = [r, 0., 0., 0., (MU / r).sqrt(), 0.];
        let (lo, hi) = orbital_bounds(y);
        assert!((lo - 500000.).abs() < 1.);
        assert!((hi - 500000.).abs() < 1.);
    }
    #[test]
    fn numerical_drag_energy_loss_and_step_convergence() {
        Atmosphere {
            altitude: vec![120., 125.],
            log_density: vec![(1e-12f64).ln(); 2],
            f107: 150.,
        }
        .install()
        .unwrap();
        let start = 1793707200.;
        let end = start + 12000.;
        let env = Environment::new(start - 1., end + 1., 150.);
        let r = RE + 400000.;
        let y = [r, 0., 0., 0., (MU / r).sqrt(), 0.];
        let a = env
            .propagate(start, end, y, 0.1, 1e-10, 120., 600.)
            .unwrap();
        let b = env.propagate(start, end, y, 0.1, 1e-11, 60., 600.).unwrap();
        let c = env.propagate(start, end, y, 0., 1e-11, 60., 600.).unwrap();
        let last = |q: &Track| -> [f64; 6] { q.state[q.state.len() - 6..].try_into().unwrap() };
        let (aa, bb, cc) = (last(&a), last(&b), last(&c));
        let energy = |z: [f64; 6]| {
            let r = norm([z[0], z[1], z[2]]);
            dot([z[3], z[4], z[5]], [z[3], z[4], z[5]]) / 2. - MU / r
        };
        assert!(energy(bb) < energy(cc) - 100.);
        let error = (0..3).map(|j| (aa[j] - bb[j]).powi(2)).sum::<f64>().sqrt();
        assert!(error < 1., "step convergence {error}m");
        let rp = RE + 500000.;
        let ra = RE + 4000000.;
        let sma = (rp + ra) / 2.;
        let high = [rp, 0., 0., 0., (MU * (2. / rp - 1. / sma)).sqrt(), 0.];
        let tr = env
            .propagate(start, end, high, 0., 1e-10, 120., 300.)
            .unwrap();
        assert!(
            tr.state
                .chunks_exact(6)
                .any(|z| norm([z[0], z[1], z[2]]) > RE + 2500000.)
        );
    }
    #[test]
    fn full_leo_grid_preserves_reentry_boundary() {
        let e = leo_altitude_edges();
        assert_eq!(e[0], 0.);
        assert_eq!(*e.last().unwrap(), 2000.);
        assert!(e.contains(&120.));
        assert!(e.windows(2).all(|w| w[1] > w[0]));
        assert!(!e.windows(2).any(|w| w[0] < 120. && w[1] > 120.));
    }
}
/// Phase mix the osculating Kepler orbit, preserving its plane and eccentricity.
pub fn mixed_state(y: [f64; 6], mean: f64) -> [f64; 6] {
    let r = [y[0], y[1], y[2]];
    let v = [y[3], y[4], y[5]];
    let rr = norm(r);
    let hh = cross(r, v);
    let h = hh.map(|x| x / norm(hh));
    let tmp = cross(v, hh);
    let ev = std::array::from_fn(|j| tmp[j] / MU - r[j] / rr);
    let ecc = norm(ev);
    let e1 = if ecc > 1e-9 {
        ev.map(|x| x / ecc)
    } else {
        r.map(|x| x / rr)
    };
    let e2 = cross(h, e1);
    let a = 1. / (2. / rr - dot(v, v) / MU);
    let mut e = mean;
    for _ in 0..15 {
        let de = (e - ecc * e.sin() - mean) / (1. - ecc * e.cos());
        e -= de;
        if de.abs() < 1e-13 {
            break;
        }
    }
    let rho = a * (1. - ecc * e.cos());
    let b = (1. - ecc * ecc).sqrt();
    let q = (MU * a).sqrt() / rho;
    let p: V = std::array::from_fn(|j| a * ((e.cos() - ecc) * e1[j] + b * e.sin() * e2[j]));
    let w: V = std::array::from_fn(|j| q * (-e.sin() * e1[j] + b * e.cos() * e2[j]));
    [p[0], p[1], p[2], w[0], w[1], w[2]]
}
