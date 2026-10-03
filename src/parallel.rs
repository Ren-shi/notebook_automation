//! Deterministic data parallelism.
//!
//! Work is split into blocks whose boundaries depend only on the problem size, never on
//! the number of threads, and partial results are combined in a fixed order. Results are
//! therefore bit-identical for any thread count and with the `parallel` feature off.

use crate::vec3::Vec3;

/// Target number of particle pairs per block in pairwise forces.
const PAIRS_PER_BLOCK: usize = 32_768;
/// Upper bound on blocks, which bounds the memory used for per-block buffers.
const MAX_BLOCKS: usize = 64;
/// Below this many particles, per-particle loops run as a single task.
#[cfg(feature = "parallel")]
const MIN_PARTICLES_PER_TASK: usize = 4_096;

/// Splits the rows `0..n` of the upper-triangular pair loop (`i < j`) into contiguous
/// ranges with roughly equal numbers of pairs. Returns a single range for small `n`.
pub(crate) fn pair_blocks(n: usize) -> Vec<(usize, usize)> {
    let total = n * n.saturating_sub(1) / 2;
    let blocks = (total / PAIRS_PER_BLOCK).clamp(1, MAX_BLOCKS);
    if blocks == 1 {
        return vec![(0, n)];
    }
    let target = total.div_ceil(blocks);
    let mut ranges = Vec::with_capacity(blocks);
    let (mut start, mut pairs) = (0, 0);
    for i in 0..n {
        pairs += n - 1 - i; // pairs in row i
        if pairs >= target {
            ranges.push((start, i + 1));
            start = i + 1;
            pairs = 0;
        }
    }
    if start < n {
        ranges.push((start, n));
    }
    ranges
}

/// Evaluates `f` on every block, in parallel when enabled; results are in block order.
pub(crate) fn map_blocks<T, F>(blocks: &[(usize, usize)], f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize, usize) -> T + Sync + Send,
{
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        blocks.par_iter().map(|&(lo, hi)| f(lo, hi)).collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        blocks.iter().map(|&(lo, hi)| f(lo, hi)).collect()
    }
}

/// Applies `f(index, item)` to every element; parallel for large slices.
pub(crate) fn for_each_indexed<T, F>(items: &mut [T], f: F)
where
    T: Send,
    F: Fn(usize, &mut T) + Sync + Send,
{
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        if items.len() >= 2 * MIN_PARTICLES_PER_TASK {
            items
                .par_iter_mut()
                .with_min_len(MIN_PARTICLES_PER_TASK)
                .enumerate()
                .for_each(|(i, x)| f(i, x));
            return;
        }
    }
    items.iter_mut().enumerate().for_each(|(i, x)| f(i, x));
}

/// Adds per-block buffers `(start, values)` into `acc`, where `values[k]` belongs to
/// particle `start + k`. Blocks are summed in order for every particle.
pub(crate) fn add_partials(acc: &mut [Vec3], partials: &[(usize, Vec<Vec3>)]) {
    for_each_indexed(acc, |i, a| {
        for (start, values) in partials {
            if i >= *start {
                *a += values[i - start];
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_blocks_cover_all_rows_with_balanced_work() {
        for n in [0, 1, 2, 100, 257, 1_000, 5_000, 20_000] {
            let blocks = pair_blocks(n);
            assert!(!blocks.is_empty() && blocks.len() <= MAX_BLOCKS);
            assert_eq!(blocks[0].0, 0);
            assert_eq!(blocks.last().unwrap().1, n);
            for w in blocks.windows(2) {
                assert_eq!(w[0].1, w[1].0);
            }
            if blocks.len() > 1 {
                let pairs = |(lo, hi): (usize, usize)| (lo..hi).map(|i| n - 1 - i).sum::<usize>();
                let max = blocks.iter().map(|&b| pairs(b)).max().unwrap();
                let total = n * (n - 1) / 2;
                // No block carries much more than its share (one extra row at most).
                assert!(max <= total / blocks.len() + n, "n = {n}");
            }
        }
    }
}
