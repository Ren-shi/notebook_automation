//! Broad-phase collision detection: candidate pairs from a uniform grid.

use std::collections::HashMap;

use crate::vec3::Vec3;

/// Pairs `(i, j)`, `i < j`, with `|x_i - x_j| < reach_i + reach_j`, in ascending order.
/// Particles with zero reach are skipped. Particles are binned into cubic cells of side
/// `2 max(reach)`, so only the 27 neighbouring cells need checking: O(N) for roughly
/// uniform sizes.
pub(crate) fn candidate_pairs(pos: &[Vec3], reach: &[f64]) -> Vec<(usize, usize)> {
    let max_reach = reach.iter().copied().fold(0.0, f64::max);
    let mut pairs = Vec::new();
    if max_reach.is_nan() || max_reach <= 0.0 {
        return pairs;
    }
    let cell = 2.0 * max_reach;
    let key = |p: Vec3| {
        (
            (p.x / cell).floor() as i64,
            (p.y / cell).floor() as i64,
            (p.z / cell).floor() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, (&p, &r)) in pos.iter().zip(reach).enumerate() {
        if r > 0.0 {
            grid.entry(key(p)).or_default().push(i);
        }
    }
    for (i, (&p, &r)) in pos.iter().zip(reach).enumerate() {
        if r <= 0.0 {
            continue;
        }
        let (cx, cy, cz) = key(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(members) = grid.get(&(cx + dx, cy + dy, cz + dz)) {
                        for &j in members {
                            if j > i && (pos[j] - p).norm() < r + reach[j] {
                                pairs.push((i, j));
                            }
                        }
                    }
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_finds_exactly_the_overlapping_pairs() {
        let mut s = 7u64;
        let mut rnd = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        let pos: Vec<Vec3> = (0..400)
            .map(|_| Vec3::new(rnd() * 10.0, rnd() * 10.0, rnd() * 10.0))
            .collect();
        let reach: Vec<f64> = (0..400)
            .map(|k| if k % 7 == 0 { 0.0 } else { 0.2 + 0.3 * rnd() })
            .collect();
        let mut brute = Vec::new();
        for i in 0..pos.len() {
            for j in i + 1..pos.len() {
                if reach[i] > 0.0
                    && reach[j] > 0.0
                    && (pos[j] - pos[i]).norm() < reach[i] + reach[j]
                {
                    brute.push((i, j));
                }
            }
        }
        assert_eq!(candidate_pairs(&pos, &reach), brute);
        assert!(!brute.is_empty());
    }
}
