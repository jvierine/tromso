//! Scientific figures for the finite-horizon debris-model experiment.
use anyhow::Result;
use hdf5::{File, types::VarLenUnicode};
use rustplotlib::*;
use serde_json::{Value, json};
use std::path::Path;
use tromso::*;
fn read(path: &str) -> Result<Vec<Value>> {
    let f = File::open(path)?;
    f.dataset("particles_json")?
        .read_raw::<VarLenUnicode>()?
        .iter()
        .map(|s| Ok(serde_json::from_str(s.as_str())?))
        .collect()
}
fn fig(h: f64) -> Result<Figure> {
    let mut f = ArticleLayout::new(174., 10.95)?.figure(1., h)?;
    f.latex.preamble =
        r"\usepackage[T1]{fontenc}\usepackage{amsmath,amssymb}\usepackage{newtxtext,newtxmath}"
            .into();
    Ok(f)
}
fn add(p: &mut Panel, pts: Vec<[f64; 2]>, j: usize, label: &str) {
    let l = Line::series(pts, j);
    let e = l.legend(label);
    p.push(Layer::Line(l));
    p.legend
        .get_or_insert_with(|| Legend::new(vec![]))
        .entries
        .push(e);
}
fn color(j: usize) -> [f64; 3] {
    let s = LineStyle::for_series(j).color;
    std::array::from_fn(|k| u8::from_str_radix(&s[1 + 2 * k..3 + 2 * k], 16).unwrap() as f64 / 255.)
}
fn main() -> Result<()> {
    std::fs::create_dir_all("output")?;
    let mut f = fig(79.)?;
    let mut ax = Axes::new(
        [0.8, 90.],
        [0.003, 2.],
        "Characteristic length / sphere diameter [mm]",
        r"Area / mass [m$^2$/kg]",
    );
    ax.set_xscale(AxisScale::Log10)?;
    ax.set_yscale(AxisScale::Log10)?;
    let mut p = Panel::new(Rect::new(25., 6., 140., 56.), ax);
    add(
        &mut p,
        (0..201)
            .map(|k| {
                let d = 0.001 * 79f64.powf(k as f64 / 200.);
                [d * 1000., sphere_am(d, 2700.)]
            })
            .collect(),
        0,
        "Compact aluminium spheres",
    );
    add(
        &mut p,
        (0..201)
            .map(|k| {
                let d = 0.001 * 79f64.powf(k as f64 / 200.);
                [d * 1000., 10f64.powf(small_am_parameters(d).unwrap().0)]
            })
            .collect(),
        1,
        "Breakup model: median small-fragment A/m",
    );
    f.panels.push(p);
    f.save_for_article(
        "output/debris-area-mass.pdf",
        "figures/debris-area-mass.pdf",
    )?;
    let mut f = fig(87.)?;
    let mut summary = vec![];
    for (col, law) in ["sphere", "sbm"].iter().enumerate() {
        let rows = read(&format!("output/lifetime-{law}-365d.h5"))?;
        let mut ax = Axes::new(
            [180., 830.],
            [0.01, 600.],
            "Initial geocentric altitude [km]",
            if col == 0 {
                "Time to 120 km, or censor [days]"
            } else {
                ""
            },
        );
        ax.set_yscale(AxisScale::Log10)?;
        ax.show_y_axis = col == 0;
        let mut panel = Panel::new(Rect::new(23. + col as f64 * 77., 7., 69., 60.), ax);
        for (j, d) in [0.001, 0.01, 0.05].iter().enumerate() {
            if law == &"sbm" && j == 1 {
                continue;
            }
            let subset: Vec<_> = rows
                .iter()
                .filter(|r| (r["characteristic_length_m"].as_f64().unwrap() - d).abs() < 1e-8)
                .collect();
            let mut good = vec![];
            let mut censored = vec![];
            for r in subset {
                let h = r["source_altitude_km"].as_f64().unwrap();
                if let Some(t) = r["reentry_after_days"].as_f64() {
                    good.push([h, t]);
                    summary.push(json!({"law":law,"size_m":d,"altitude_km":h,"lifetime_days":t,"censored":false}));
                } else {
                    censored.push([h, 365.]);
                    summary.push(json!({"law":law,"size_m":d,"altitude_km":h,"lifetime_days_lower_bound":365.,"censored":true}));
                }
            }
            let label = if law == &"sbm" && j == 0 {
                "1 mm / 1 cm"
            } else if j == 0 {
                "1 mm"
            } else if j == 1 {
                "1 cm"
            } else {
                "5 cm"
            };
            add(&mut panel, good.clone(), j, label);
            panel.markers(
                &mut f,
                &good,
                &MarkerStyle {
                    shape: MarkerShape::Circle,
                    size_pt: 4.,
                    color: color(j),
                },
            )?;
            panel.markers(
                &mut f,
                &censored,
                &MarkerStyle {
                    shape: MarkerShape::Triangle,
                    size_pt: 5.,
                    color: color(j),
                },
            )?;
        }
        f.panels.push(panel);
    }
    f.save_for_article(
        "output/debris-lifetimes.pdf",
        "figures/debris-lifetimes.pdf",
    )?;
    let out = File::create("output/plot-summary.h5")?;
    out.new_dataset_builder()
        .with_data(serde_json::to_string(&summary)?.as_bytes())
        .create("lifetime_summary_json_utf8")?;
    println!("{}", serde_json::to_string_pretty(&summary)?);
    let pop = "output/catalogue-sbm-365d.h5";
    if Path::new(pop).exists() {
        let rows = read(pop)?;
        let mut f = fig(81.)?;
        let mut ax = Axes::new(
            [100., 2100.],
            [1e-14, 1e-8],
            "Current geocentric altitude [km]",
            r"Weighted objects / km$^3$ / $\ln L_c$",
        );
        ax.set_yscale(AxisScale::Log10)?;
        let mut p = Panel::new(Rect::new(26., 6., 139., 57.), ax);
        let groups = [
            (0.001, 0.003, "1--3 mm"),
            (0.003, 0.01, "3--10 mm"),
            (0.01, 0.079, "1--7.9 cm"),
        ];
        let mut density_rows = vec![];
        for (j, (lo, hi, label)) in groups.iter().enumerate() {
            let mut pts = vec![];
            for k in 0..10 {
                let h0 = 100. + k as f64 * 200.;
                let h1 = h0 + 200.;
                let vol = 4. * std::f64::consts::PI / 3.
                    * (((RE / 1000. + h1).powi(3)) - (RE / 1000. + h0).powi(3));
                let w: f64 = rows
                    .iter()
                    .filter(|r| {
                        r["status"] == "surviving_censored"
                            && r["characteristic_length_m"].as_f64().unwrap() >= *lo
                            && r["characteristic_length_m"].as_f64().unwrap() < *hi
                            && r["final_altitude_km"].as_f64().unwrap() >= h0
                            && r["final_altitude_km"].as_f64().unwrap() < h1
                    })
                    .map(|r| r["weight_objects"].as_f64().unwrap())
                    .sum();
                let z = w / (vol * (hi / lo).ln());
                density_rows.push(json!({"size_min_m":lo,"size_max_m":hi,"altitude_min_km":h0,"altitude_max_km":h1,"density_per_km3_per_ln_size":z}));
                if z > 0. {
                    pts.push([0.5 * (h0 + h1), z]);
                } else {
                    pts.push([0.5 * (h0 + h1), f64::NAN]);
                }
            }
            add(&mut p, pts, j, label);
        }
        f.panels.push(p);
        f.save_for_article(
            "output/debris-size-altitude.pdf",
            "figures/debris-size-altitude.pdf",
        )?;
        out.new_dataset_builder()
            .with_data(serde_json::to_string(&density_rows)?.as_bytes())
            .create("density_summary_json_utf8")?;
    }
    Ok(())
}
