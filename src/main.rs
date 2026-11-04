//! Weighted source-to-survival calculation using the owner's DOP853 propagator.
use anyhow::{Result, ensure};
use clap::Parser;
use hdf5::types::VarLenUnicode;
use rayon::prelude::*;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{path::PathBuf, time::Instant};
use tromso::*;
#[derive(Parser, Debug)]
struct Args {
    #[arg(long, default_value = "2026-11-03T12:00:00Z")]
    epoch: String,
    #[arg(long)]
    catalogue: PathBuf,
    #[arg(long, default_value = "output/population.h5")]
    output: PathBuf,
    #[arg(long, default_value_t = 384)]
    particles: usize,
    #[arg(long, default_value_t = 30.)]
    history_days: f64,
    #[arg(long, default_value = "sbm-small")]
    area_mass: String,
    #[arg(long, default_value = "collision")]
    event: String,
    #[arg(long, default_value = "all")]
    parents: String,
    #[arg(long, default_value_t = 150.)]
    f107: f64,
    #[arg(long, default_value_t = 2.2)]
    cd: f64,
    #[arg(long, default_value_t = 1.)]
    births_per_day: f64,
    #[arg(long, default_value_t = 6)]
    threads: usize,
    #[arg(long, default_value_t = 20261104)]
    seed: u64,
    #[arg(long, default_value_t = 1e-9)]
    rtol: f64,
    #[arg(long, default_value_t = 120.)]
    max_step_s: f64,
    #[arg(long, default_value_t = 1800.)]
    sample_step_s: f64,
    #[arg(long, default_value_t = false)]
    lifetime_grid: bool,
}
fn number(v: &Value, k: &str) -> Result<f64> {
    v[k].as_f64()
        .or_else(|| v[k].as_str().and_then(|s| s.parse().ok()))
        .ok_or_else(|| anyhow::anyhow!("Missing numeric {k}"))
}
fn strv<'a>(v: &'a Value, k: &str) -> &'a str {
    v[k].as_str().unwrap_or("")
}
fn parents(a: &Args) -> Result<(Vec<Parent>, Value)> {
    let f = hdf5::File::open(&a.catalogue)?;
    let rows = f.dataset("records_json")?.read_raw::<VarLenUnicode>()?;
    let total = rows.len();
    let mut out = vec![];
    let mut failures = 0;
    for s in rows {
        let mut row: Value = serde_json::from_str(s.as_str())?;
        let typ = strv(&row, "OBJECT_TYPE").to_owned();
        if a.parents == "debris" && typ != "DEBRIS" {
            continue;
        }
        let mm = number(&row, "MEAN_MOTION")?;
        let ecc = number(&row, "ECCENTRICITY")?;
        if !(mm > 0. && ecc < 1.) {
            continue;
        }
        let sma = (MU / (mm * std::f64::consts::TAU / DAY).powi(2)).cbrt();
        if sma * (1. - ecc) - RE < 120000. || sma * (1. - ecc) - RE > 2000000. {
            continue;
        }
        let t = chrono::NaiveDateTime::parse_from_str(
            strv(&row, "EPOCH").trim_end_matches('Z'),
            "%Y-%m-%dT%H:%M:%S%.f",
        )?
        .and_utc()
        .timestamp() as f64;
        for (k, v) in [
            ("CLASSIFICATION_TYPE", json!("U")),
            ("OBJECT_ID", json!("")),
            ("EPHEMERIS_TYPE", json!(0)),
            ("ELEMENT_SET_NO", json!(0)),
            ("REV_AT_EPOCH", json!(0)),
            ("MEAN_MOTION_DOT", json!(0)),
            ("MEAN_MOTION_DDOT", json!(0)),
        ] {
            if row[k].is_null() {
                row[k] = v;
            }
        }
        row["EPOCH"] = json!(strv(&row, "EPOCH").trim_end_matches('Z'));
        let el: sgp4::Elements = serde_json::from_value(row.clone())?;
        let c = sgp4::Constants::new(
            sgp4::WGS72,
            sgp4::iau_epoch_to_sidereal_time,
            el.epoch(),
            el.drag_term,
            sgp4::Orbit::from_kozai_elements(
                &sgp4::WGS72,
                el.inclination.to_radians(),
                el.right_ascension.to_radians(),
                el.eccentricity,
                el.argument_of_perigee.to_radians(),
                el.mean_anomaly.to_radians(),
                el.mean_motion * std::f64::consts::PI / 720.,
            )?,
        )?;
        match c.propagate(sgp4::MinutesSinceEpoch(0.)) {
            Ok(p) => {
                let y = teme_gcrs(t, p.position, p.velocity);
                let (lo, _) = orbital_bounds(y);
                if lo > 120000. && lo < 2100000. {
                    out.push(Parent {
                        id: number(&row, "NORAD_CAT_ID")? as u64,
                        name: strv(&row, "OBJECT_NAME").into(),
                        epoch: t,
                        state: y,
                        object_type: typ,
                    });
                }
            }
            Err(_) => {
                failures += 1;
            }
        }
    }
    ensure!(!out.is_empty(), "No eligible LEO catalogue parents");
    let n = out.len();
    Ok((
        out,
        json!({"catalogue_rows":total,"selected_parents":n,"sgp4_failures":failures,"selection":"mean-element perigee120--2000km; osculating perigee120--2100km; no apogee ceiling, includes LEO-crossing eccentric orbits","parent_filter":a.parents,"birth_kernel":"uniform selected catalogue objects, phase-mixed osculating orbit; stationary source geometry reused at each birth time, NOT reconstructed event history","completeness_correction":false}),
    ))
}
#[derive(serde::Serialize)]
struct Row {
    parent: u64,
    size_bin: usize,
    characteristic_length_m: f64,
    area_mass_m2_kg: f64,
    mass_kg: f64,
    age_days: f64,
    source_altitude_km: f64,
    weight_objects: f64,
    delta_v_m_s: f64,
    initial: [f64; 6],
    final_state: [f64; 6],
    final_altitude_km: f64,
    perigee_km: f64,
    status: String,
    reentry_after_days: Option<f64>,
    nfev: usize,
    error: Option<String>,
}
fn save(path: &PathBuf, rows: &[Row], meta: &Value, atm: &Atmosphere) -> Result<()> {
    if let Some(p) = path.parent() {
        std::fs::create_dir_all(p)?;
    }
    let tmp = path.with_extension("partial.h5");
    let f = hdf5::File::create(&tmp)?;
    let strings: Vec<_> = rows
        .iter()
        .map(|r| {
            serde_json::to_string(r)
                .unwrap()
                .parse::<VarLenUnicode>()
                .unwrap()
        })
        .collect();
    f.new_dataset_builder()
        .with_data(&strings)
        .create("particles_json")?;
    f.new_dataset_builder()
        .with_data(meta.to_string().as_bytes())
        .create("metadata_json_utf8")?;
    let g = f.create_group("atmosphere")?;
    g.new_dataset_builder()
        .with_data(&atm.altitude)
        .create("altitude_km")?;
    g.new_dataset_builder()
        .with_data(&atm.log_density)
        .create("log_density_kg_m3")?;
    let sizes: Vec<_> = rows.iter().map(|r| r.characteristic_length_m).collect();
    f.new_dataset_builder()
        .with_data(&sizes)
        .create("characteristic_length_m")?;
    let weights: Vec<_> = rows.iter().map(|r| r.weight_objects).collect();
    f.new_dataset_builder()
        .with_data(&weights)
        .create("population_weight")?;
    let states: Vec<_> = rows.iter().flat_map(|r| r.final_state).collect();
    f.new_dataset::<f64>()
        .shape((rows.len(), 6))
        .create("final_state_gcrs")?
        .write_raw(&states)?;
    let edges = leo_altitude_edges();
    let sizes: Vec<_> = (0..=6).map(|j| 0.001 * 79f64.powf(j as f64 / 6.)).collect();
    let nh = edges.len() - 1;
    let ns = 6;
    let mut population = vec![0.; nh * ns];
    let mut w2 = vec![0.; nh * ns];
    let mut counts = vec![0u64; nh * ns];
    for row in rows.iter().filter(|r| r.status == "surviving_censored") {
        let h = row.final_altitude_km;
        let d = row.characteristic_length_m;
        if h >= 120. && h < 2000. && d >= sizes[0] && d < sizes[ns] {
            let ih = edges.partition_point(|x| *x <= h) - 1;
            let is = sizes.partition_point(|x| *x <= d) - 1;
            let k = ih * ns + is;
            population[k] += row.weight_objects;
            w2[k] += row.weight_objects.powi(2);
            counts[k] += 1;
        }
    }
    let valid: Vec<u8> = edges[..nh].iter().map(|h| (*h >= 120.) as u8).collect();
    let mut density = population.clone();
    for ih in 0..nh {
        let volume = 4. * std::f64::consts::PI / 3.
            * ((RE / 1000. + edges[ih + 1]).powi(3) - (RE / 1000. + edges[ih]).powi(3));
        for is in 0..ns {
            let k = ih * ns + is;
            if valid[ih] == 0 {
                population[k] = f64::NAN;
                w2[k] = f64::NAN;
                density[k] = f64::NAN;
            } else {
                density[k] /= volume * (sizes[is + 1] / sizes[is]).ln();
            }
        }
    }
    let g = f.create_group("leo")?;
    g.new_dataset_builder()
        .with_data(&edges)
        .create("altitude_edges_km")?;
    g.new_dataset_builder()
        .with_data(&sizes)
        .create("characteristic_length_edges_m")?;
    g.new_dataset_builder()
        .with_data(&valid)
        .create("orbital_population_valid")?;
    for (name, data) in [
        ("weighted_objects", &population),
        ("weighted_sum_squares", &w2),
        ("density_per_km3_per_ln_length", &density),
    ] {
        g.new_dataset::<f64>()
            .shape((nh, ns))
            .create(name)?
            .write_raw(data)?;
    }
    g.new_dataset::<u64>()
        .shape((nh, ns))
        .create("monte_carlo_counts")?
        .write_raw(&counts)?;
    f.flush()?;
    drop(f);
    std::fs::rename(tmp, path)?;
    Ok(())
}
fn main() -> Result<()> {
    let a = Args::parse();
    ensure!(
        a.particles >= 6
            && a.history_days > 0.
            && a.cd > 0.
            && a.f107 > 0.
            && a.births_per_day >= 0.,
        "Invalid model parameters"
    );
    ensure!(
        ["collision", "explosion"].contains(&a.event.as_str())
            && ["sphere", "sbm-small"].contains(&a.area_mass.as_str())
            && ["all", "debris"].contains(&a.parents.as_str()),
        "Model selector"
    );
    rayon::ThreadPoolBuilder::new()
        .num_threads(a.threads)
        .build_global()?;
    let timer = Instant::now();
    let (ps, parent_info) = parents(&a)?;
    eprintln!(
        "{} catalogue parents; building atmospheric climatology",
        ps.len()
    );
    let atm = Atmosphere::build(a.f107);
    atm.clone().install()?;
    let now = chrono::DateTime::parse_from_rfc3339(&a.epoch)?.timestamp() as f64;
    let env = Environment::new(now - a.history_days * DAY - DAY, now + DAY, a.f107);
    let p = if a.event == "collision" { 1.71 } else { 1.6 };
    let edges: Vec<_> = (0..=6)
        .map(|j| 0.001 * (79f64).powf(j as f64 / 6.))
        .collect();
    let nb = 6;
    let mut jobs = vec![];
    if a.lifetime_grid {
        for height in [150., 200., 300., 400., 600., 800., 1000., 1500., 2000.] {
            for d in [0.001, 0.01, 0.05] {
                jobs.push((height, d));
            }
        }
    }
    let count = if a.lifetime_grid {
        jobs.len()
    } else {
        a.particles
    };
    let completed = std::sync::atomic::AtomicUsize::new(0);
    let rows: Vec<Row> = (0..count)
        .into_par_iter()
        .map(|k| {
            let mut rng = Rng(a.seed.wrapping_add(k as u64 * 0x100000001b3));
            let bin = k % nb;
            let mut q = Row {
                parent: 0,
                size_bin: bin,
                characteristic_length_m: 0.,
                area_mass_m2_kg: 0.,
                mass_kg: 0.,
                age_days: 0.,
                source_altitude_km: 0.,
                weight_objects: 0.,
                delta_v_m_s: 0.,
                initial: [0.; 6],
                final_state: [0.; 6],
                final_altitude_km: 0.,
                perigee_km: 0.,
                status: "failed".into(),
                reentry_after_days: None,
                nfev: 0,
                error: None,
            };
            let result = (|| -> Result<()> {
                let (y, d, age, w) = if a.lifetime_grid {
                    let (h, d) = jobs[k];
                    let r = RE + h * 1000.;
                    let inc = 53f64.to_radians();
                    (
                        [
                            r,
                            0.,
                            0.,
                            0.,
                            (MU / r).sqrt() * inc.cos(),
                            (MU / r).sqrt() * inc.sin(),
                        ],
                        d,
                        a.history_days,
                        1.,
                    )
                } else {
                    let parent = &ps[(rng.uniform() * ps.len() as f64) as usize];
                    q.parent = parent.id;
                    let y = mixed_state(parent.state, std::f64::consts::TAU * rng.uniform());
                    let d = power_quantile(rng.uniform(), edges[bin], edges[bin + 1], p);
                    let age = a.history_days * rng.uniform();
                    let nbin = (0..count).filter(|j| j % nb == bin).count();
                    let w = a.births_per_day
                        * a.history_days
                        * bin_probability(edges[bin], edges[bin + 1], edges[0], edges[nb], p)
                        / nbin as f64;
                    (y, d, age, w)
                };
                q.characteristic_length_m = d;
                q.age_days = age;
                q.weight_objects = w;
                q.area_mass_m2_kg = if a.lifetime_grid && a.area_mass == "sbm-small" {
                    10f64.powf(small_am_parameters(d)?.0)
                } else {
                    sample_am(d, &a.area_mass, &mut rng)?
                };
                q.mass_kg = if a.area_mass == "sphere" {
                    2700. * std::f64::consts::PI * d.powi(3) / 6.
                } else {
                    area(d) / q.area_mass_m2_kg
                };
                q.source_altitude_km = (norm([y[0], y[1], y[2]]) - RE) / 1000.;
                let dv = if a.lifetime_grid {
                    [0.; 3]
                } else {
                    ejection(q.area_mass_m2_kg, a.event == "collision", &mut rng)
                };
                q.delta_v_m_s = norm(dv);
                let mut initial = y;
                for j in 0..3 {
                    initial[j + 3] += dv[j];
                }
                q.initial = initial;
                q.final_state = initial;
                let (peri, apo) = orbital_bounds(initial);
                q.perigee_km = peri / 1000.;
                if !apo.is_finite() {
                    q.status = "unbound_ejection".into();
                    return Ok(());
                }
                if norm([initial[0], initial[1], initial[2]]) < RE + 120000. {
                    q.status = "born_below_stop".into();
                    return Ok(());
                }
                let birth = now - age * DAY;
                let tr = env.propagate(
                    birth,
                    now,
                    initial,
                    a.cd * q.area_mass_m2_kg,
                    a.rtol,
                    a.max_step_s,
                    a.sample_step_s,
                )?;
                q.nfev = tr.nfev;
                q.final_state = tr.state[tr.state.len() - 6..].try_into().unwrap();
                q.final_altitude_km =
                    (norm([q.final_state[0], q.final_state[1], q.final_state[2]]) - RE) / 1000.;
                if tr.event.is_finite() {
                    q.status = "reentered_120km".into();
                    q.reentry_after_days = Some((tr.event - birth) / DAY);
                } else {
                    q.status = "surviving_censored".into();
                }
                Ok(())
            })();
            if let Err(e) = result {
                q.error = Some(e.to_string());
            }
            let done = completed.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
            if done % 32 == 0 || done == count {
                eprintln!("{done}/{count} trajectories complete");
            }
            q
        })
        .collect();
    let failures = rows.iter().filter(|r| r.status == "failed").count();
    let survived = rows
        .iter()
        .filter(|r| r.status == "surviving_censored")
        .count();
    let reentered = rows
        .iter()
        .filter(|r| r.reentry_after_days.is_some())
        .count();
    let weighted: f64 = rows
        .iter()
        .filter(|r| r.status == "surviving_censored")
        .map(|r| r.weight_objects)
        .sum();
    let mut meta = json!({"schema":"tromso_source_survival_v2","parent_info":parent_info,"catalogue_sha256":format!("{:x}",Sha256::digest(std::fs::read(&a.catalogue)?)),"snapshot_geometry_reference":"2026-11-03","evaluation_epoch_utc":"2026-11-03T12:00:00Z","particles":count,"history_days":a.history_days,"births_per_day_1mm_to79mm":a.births_per_day,"weighted_surviving_objects":weighted,"survived":survived,"reentered":reentered,"failed":failures,"cumulative_size_exponent":p,"size_min_m":0.001,"size_max_m":0.079,"size_coordinate":"NASA characteristic length Lc; sphere sensitivity case identifies it with diameter","area_mass_model":a.area_mass,"ejection":"isotropic direction; log10 speed N(0.9 log10 A/m+2.9,0.4) collision or N(0.2 log10 A/m+1.85,0.4) explosion","ejection_truncation":"none; unbound and immediate-loss trajectories retained explicitly","birth_age_prior":"uniform over finite source history, stationary catalogue geometry; no reconstruction of historical collisions","normalization":"specified source rate, NOT fitted absolute debris population; does not use a catalogue count as number of breakup events","forces":"owner numerical_propagator DOP853, central gravity+J2+SOFA Sun/Moon third bodies, rigidly corotating atmosphere, Cd A/m; no SRP or fragment-fragment collisions","density":"pure Rust NRLMSISE-00 gtd7d arithmetic mean over 5 uniform-sin-lat nodes,4UT,2equinox days; log-height table every5km120-2500km; constant F107,Ap4; not observed weather history","f107_sfu":a.f107,"cd":a.cd,"frame":"GCRS initial state from TEME with IAU2006/2000A, TT-UTC69.184s, UT1=UTC, polar motion0; not precise EOP","numerical_rtol":a.rtol,"max_step_s":a.max_step_s,"stop_geodetic_km":120.,"lifetime_grid":a.lifetime_grid,"runtime_s":timer.elapsed().as_secs_f64(),"published_model_scope":"small-fragment A/m prescription restricted to Lc<8cm; not full NASA breakup event mass conservation or calibrated ORDEM/MASTER environment","sources":["https://doi.org/10.1007/s10569-024-10205-3","https://doi.org/10.1016/S0273-1177(01)00423-9","https://doi.org/10.1029/2002JA009430"]});
    let mut status_weights = std::collections::BTreeMap::<String, f64>::new();
    for row in &rows {
        *status_weights.entry(row.status.clone()).or_default() += row.weight_objects;
    }
    meta["weighted_status_objects"] = json!(status_weights);
    meta["random_seed"] = json!(a.seed);
    meta["evaluation_epoch_utc"] = json!(&a.epoch);
    meta["snapshot_geometry_reference"] =
        json!("individual GP epochs retained; see catalogue hash and source-parent epoch range");
    meta["source_parent_epoch_min_unix"] =
        json!(ps.iter().map(|p| p.epoch).fold(f64::INFINITY, f64::min));
    meta["source_parent_epoch_max_unix"] =
        json!(ps.iter().map(|p| p.epoch).fold(f64::NEG_INFINITY, f64::max));
    meta["threads"] = json!(a.threads);
    meta["sample_step_s"] = json!(a.sample_step_s);
    meta["drag_momentum_coefficient_gamma"] = json!(a.cd / 2.);
    meta["reporting_domain_km"] = json!([0, 2000]);
    meta["below_120km"] =
        json!("reentry/ablation regime not propagated: masked NaN, not an inferred zero density");
    meta["high_apogee_fragments"] = json!(
        "bound fragments retained without apogee exclusion; external motion propagated and only residence in0--2000km reported"
    );
    meta["upper_atmosphere"] = json!("drag neglected above2500km; trajectories not removed there");
    save(&a.output, &rows, &meta, &atm)?;
    println!("{}", serde_json::to_string_pretty(&meta)?);
    ensure!(
        failures == 0,
        "{failures} failed trajectories retained in output"
    );
    Ok(())
}
