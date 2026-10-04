//! The nuclear event generator is reproducible: the same seed gives bit-identical events for
//! any thread count and however the run is split into pieces.

use physim::nuclear::{Generator, Record};

fn bits(records: &[Record]) -> Vec<u64> {
    records
        .iter()
        .flat_map(|r| {
            [
                r.event,
                r.detector as u64,
                r.segment.0 as u64 * 1000 + r.segment.1 as u64,
                r.recoil as u64,
                r.measured.to_bits(),
                r.theta.to_bits(),
                r.phi.to_bits(),
                r.depth.to_bits(),
                r.weight.to_bits(),
            ]
        })
        .collect()
}

#[test]
fn events_are_identical_for_any_thread_count() {
    let gen = Generator::demo();
    let mut results = Vec::new();
    for threads in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        results.push(pool.install(|| bits(&gen.run(50_000, 50_000, 7, 0).unwrap())));
    }
    assert!(!results[0].is_empty());
    assert_eq!(results[0], results[1]);
    let other = bits(&gen.run(50_000, 50_000, 8, 0).unwrap());
    assert_ne!(other, results[0]);
}

#[test]
fn a_run_split_in_pieces_gives_the_same_events() {
    let gen = Generator::demo();
    let whole = gen.run(30_000, 30_000, 3, 0).unwrap();
    let mut parts = gen.run(10_000, 30_000, 3, 0).unwrap();
    parts.extend(gen.run(20_000, 30_000, 3, 10_000).unwrap());
    assert_eq!(bits(&whole), bits(&parts));
}

#[test]
fn weights_add_up_to_the_sampled_rate() {
    // Every event's weight stands for a share of the channel's rate into its sampled angle
    // range; the ejectiles in the large backward annulus and the discs see a fraction of it.
    let gen = Generator::demo();
    let n = 200_000;
    let records = gen.run(n, n, 11, 0).unwrap();
    assert!(records
        .iter()
        .all(|r| r.weight > 0.0 && r.weight.is_finite()));
    assert!(records.iter().any(|r| r.detector == 3 && r.segment.0 == 15));
    // Rate at 5.5 MeV (CM 5.39 MeV) into 5°-87.5° is ~2.5e5/s; the detectors see a few percent.
    let seen: f64 = records.iter().filter(|r| !r.recoil).map(|r| r.weight).sum();
    assert!(seen > 1e3 && seen < 2.5e5, "{seen}");
    let mut bad = gen.clone();
    bad.channels[0].u_min = 0.0;
    assert!(bad.run(10, 10, 1, 0).is_err());
}

#[test]
fn excitation_lowers_energies_and_its_probability_weights_events() {
    let elastic = Generator::demo();
    let mut excited = elastic.clone();
    excited.channels[0].excitation = 0.5;
    excited.channels[0].excite_recoil = true;
    excited.channels[0].p_table = vec![0.5; 181];
    let (a, b) = (
        elastic.run(20_000, 20_000, 2, 0).unwrap(),
        excited.run(20_000, 20_000, 2, 0).unwrap(),
    );
    // Same draws, so the same events: each weight is halved, and the scattered α particles leave
    // with less energy. (Lab angles shift slightly too, so compare events seen in both runs.)
    let seen: std::collections::HashMap<u64, &physim::nuclear::Record> = a
        .iter()
        .filter(|r| !r.recoil)
        .map(|r| (r.event, r))
        .collect();
    let mut compared = 0;
    for y in b.iter().filter(|r| !r.recoil) {
        if let Some(x) = seen.get(&y.event) {
            assert!(
                (y.weight / x.weight - 0.5).abs() < 1e-3,
                "{} {}",
                x.weight,
                y.weight
            );
            assert!(y.energy < x.energy - 0.4);
            compared += 1;
        }
    }
    assert!(compared > 100, "{compared}");
    let mut bad = excited.clone();
    bad.channels[0].p_table[3] = -1.0;
    assert!(bad.run(10, 10, 1, 0).is_err());
}
