//! Reduced run audit; explicitly retain source weight outside LEO force scope.
use anyhow::Result;
use hdf5::{File, types::VarLenUnicode};
use serde_json::{Value, json};
use std::collections::BTreeMap;
fn main() -> Result<()> {
    let mut all = vec![];
    for name in [
        "catalogue-sbm-365d",
        "catalogue-sphere-365d",
        "lifetime-sphere-365d",
        "lifetime-sbm-365d",
    ] {
        let path = format!("output/{name}.h5");
        if !std::path::Path::new(&path).exists() {
            continue;
        }
        let f = File::open(path)?;
        let rows: Vec<Value> = f
            .dataset("particles_json")?
            .read_raw::<VarLenUnicode>()?
            .iter()
            .map(|s| serde_json::from_str(s.as_str()).unwrap())
            .collect();
        let meta: Value =
            serde_json::from_slice(&f.dataset("metadata_json_utf8")?.read_raw::<u8>()?)?;
        let mut counts = BTreeMap::<String, u64>::new();
        let mut weights = BTreeMap::<String, f64>::new();
        for r in &rows {
            let status = r["status"].as_str().unwrap().to_owned();
            *counts.entry(status.clone()).or_default() += 1;
            *weights.entry(status).or_default() += r["weight_objects"].as_f64().unwrap();
        }
        all.push(json!({"run":name,"counts":counts,"weighted_objects":weights,"source_history_days":meta["history_days"],"catalogue_sha256":meta["catalogue_sha256"],"parent_info":meta["parent_info"]}));
    }
    let f = File::create("output/model-summary.h5")?;
    f.new_dataset_builder()
        .with_data(serde_json::to_string(&all)?.as_bytes())
        .create("summary_json_utf8")?;
    println!("{}", serde_json::to_string_pretty(&all)?);
    Ok(())
}
