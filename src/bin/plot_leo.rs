//! Full-domain cell map, preserving invalid/reentry cells and unsampled bins.
use anyhow::{Result, ensure};
use rustplotlib::*;
fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .unwrap_or("output/full-leo-30d.h5".into());
    let f = hdf5::File::open(path)?;
    let g = f.group("leo")?;
    let x = g
        .dataset("characteristic_length_edges_m")?
        .read_raw::<f64>()?
        .into_iter()
        .map(|d| 1000. * d)
        .collect::<Vec<_>>();
    let y = g.dataset("altitude_edges_km")?.read_raw::<f64>()?;
    ensure!(
        y[0] == 0. && *y.last().unwrap() == 2000.,
        "Full LEO domain required"
    );
    let values = g
        .dataset("density_per_km3_per_ln_length")?
        .read_raw::<f64>()?
        .into_iter()
        .map(|v| {
            if v.is_finite() && v > 0. {
                v.log10()
            } else {
                f64::NAN
            }
        })
        .collect::<Vec<_>>();
    let finite = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();
    ensure!(!finite.is_empty(), "No populated LEO cells");
    let lo = finite.iter().copied().fold(f64::INFINITY, f64::min).floor();
    let hi = finite
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
        .ceil();
    let mut fig = ArticleLayout::new(174., 10.95)?.figure(1., 108.)?;
    fig.latex.preamble =
        r"\usepackage[T1]{fontenc}\usepackage{amsmath,amssymb}\usepackage{newtxtext,newtxmath}"
            .into();
    let s = fig.scale(Scale::new(
        [lo, hi.max(lo + 1.)],
        Normalization::Linear,
        ColorMap::sequential(),
    )?);
    let mut ax = Axes::new(
        [x[0], *x.last().unwrap()],
        [0., 2000.],
        "Characteristic length [mm]",
        "Geocentric altitude [km]",
    );
    ax.set_xscale(AxisScale::Log10)?;
    let mut p = Panel::new(Rect::new(22., 6., 119., 84.), ax);
    p.push(Layer::Pcolormesh(Field::new(x.clone(), y, values, s)?));
    let mut line = Line::new(vec![[x[0], 120.], [*x.last().unwrap(), 120.]], "#444444");
    line.dash = Some("1.2 1.2".into());
    p.push(Layer::Line(line));
    fig.panels.push(p);
    fig.colorbars.push(Colorbar::new(
        Rect::new(147., 6., 3., 84.),
        s,
        r"$\log_{10}$ weighted density [km$^{-3}$ / $\ln L_c$]",
        (lo as i32..=hi as i32).map(|i| i as f64).collect(),
    ));
    std::fs::create_dir_all("output")?;
    fig.save_for_article("output/tromso-full-leo.pdf", "figures/tromso-full-leo.pdf")?;
    Ok(())
}
