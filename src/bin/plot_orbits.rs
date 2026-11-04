//! Orthographic catalogue view. Point diameters are illustrative, not object sizes.
use anyhow::Result;
use hdf5::types::VarLenUnicode;
use rustplotlib::*;
use serde_json::{Value, json};
fn project(p: [f64; 3]) -> ([f64; 2], f64) {
    let a = 0.55_f64;
    let b = 0.40_f64;
    let x = a.cos() * p[0] + a.sin() * p[1];
    let y = -a.sin() * p[0] + a.cos() * p[1];
    (
        [x, b.sin() * y + b.cos() * p[2]],
        b.cos() * y - b.sin() * p[2],
    )
}
fn main() -> Result<()> {
    let input = std::env::args()
        .nth(1)
        .unwrap_or("data/catalogue-2026-11-03.h5".into());
    let file = hdf5::File::open(input)?;
    let epoch = chrono::DateTime::parse_from_rfc3339("2026-11-03T12:00:00Z")?.timestamp() as f64;
    let mut points = Vec::new();
    let mut failures = 0;
    for s in file.dataset("records_json")?.read_raw::<VarLenUnicode>()? {
        let mut row: Value = serde_json::from_str(s.as_str())?;
        let t = chrono::NaiveDateTime::parse_from_str(
            row["EPOCH"].as_str().unwrap().trim_end_matches('Z'),
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
        row["EPOCH"] = json!(row["EPOCH"].as_str().unwrap().trim_end_matches('Z'));
        let result = (|| -> Result<[f64; 3]> {
            let e: sgp4::Elements = serde_json::from_value(row)?;
            let c = sgp4::Constants::from_elements(&e)?;
            let p = c.propagate(sgp4::MinutesSinceEpoch((epoch - t) / 60.))?;
            let y = tromso::teme_gcrs(epoch, p.position, p.velocity);
            Ok([y[0] / 1000., y[1] / 1000., y[2] / 1000.])
        })();
        match result {
            Ok(p) => points.push(p),
            Err(_) => failures += 1,
        }
    }
    let mut fig = ArticleLayout::new(174., 10.95)?.figure(1., 104.)?;
    fig.latex.preamble = r"\usepackage[T1]{fontenc}\usepackage{newtxtext,newtxmath}".into();
    for (i, lim) in [51000., 10500.].into_iter().enumerate() {
        let mut panel = Panel::new(
            Rect::new(3. + i as f64 * 86., 8., 82., 82.),
            Axes::new([-lim, lim], [-lim, lim], "", ""),
        );
        panel.axes_visible = false;
        // Draw far-side objects, opaque Earth, then foreground objects.
        for front in [false, true] {
            let mut pts = vec![];
            for &p in &points {
                let (q, z) = project(p);
                if q[0].abs() > lim || q[1].abs() > lim {
                    continue;
                }
                if (z >= 0.) == front {
                    pts.push(q);
                }
            }
            panel.markers(
                &mut fig,
                &pts,
                &MarkerStyle {
                    shape: MarkerShape::Circle,
                    size_pt: if i == 0 { 0.48 } else { 0.65 },
                    color: [0., 0.447, 0.698],
                },
            )?;
            if !front {
                panel.markers(
                    &mut fig,
                    &[[0., 0.]],
                    &MarkerStyle {
                        shape: MarkerShape::Circle,
                        size_pt: 2. * 6378.137 / (2. * lim) * 82. * 72.27 / 25.4,
                        color: [0.79, 0.85, 0.88],
                    },
                )?;
                for lat in [-60_f64, -30., 0., 30., 60.] {
                    let lat = lat.to_radians();
                    let mut seg = vec![];
                    for k in 0..=360 {
                        let l = (k as f64).to_radians();
                        let (q, z) = project([
                            6378.137 * lat.cos() * l.cos(),
                            6378.137 * lat.cos() * l.sin(),
                            6378.137 * lat.sin(),
                        ]);
                        if z >= 0. {
                            seg.push(q)
                        } else if !seg.is_empty() {
                            panel.push(Layer::Line(Line::new(std::mem::take(&mut seg), "#ffffff")));
                        }
                    }
                    if !seg.is_empty() {
                        panel.push(Layer::Line(Line::new(seg, "#ffffff")));
                    }
                }
            }
        }
        panel.label_at(
            if i == 0 {
                "LEO to geostationary orbit"
            } else {
                "Low Earth orbit"
            },
            [0., lim * 1.13],
        );
        fig.panels.push(panel);
    }
    std::fs::create_dir_all("output")?;
    fig.save_for_article("output/earth-debris.pdf", "figures/earth-debris.pdf")?;
    println!(
        "{} propagated catalogue objects; {} failures",
        points.len(),
        failures
    );
    Ok(())
}
