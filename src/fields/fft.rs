//! A small radix-2 complex FFT and its multi-dimensional form on row-major grids.
//!
//! Along the leading axes the butterflies act on whole contiguous rows of the last axis at
//! once (cache-friendly and vectorisable); along the last axis they act on single lines.
//! Lines can be restricted to a box of indices ([`fft_nd_pruned`]), which lets zero-padded
//! convolutions skip the lines that are known to be zero or are not needed.

use crate::error::{invalid, Result};

/// A complex number as `(re, im)`.
pub type Complex = (f64, f64);

/// `e^{-2πik/n}` for `k < n/2`.
fn twiddles(n: usize) -> Vec<Complex> {
    (0..n / 2)
        .map(|k| {
            let (s, c) = (-2.0 * std::f64::consts::PI * k as f64 / n as f64).sin_cos();
            (c, s)
        })
        .collect()
}

fn bit_reverse_index(n: usize) -> Vec<usize> {
    let bits = n.trailing_zeros();
    (0..n)
        .map(|i| {
            if bits == 0 {
                0
            } else {
                i.reverse_bits() >> (usize::BITS - bits)
            }
        })
        .collect()
}

/// FFT along one axis of rows: element `e` of position `t` is at `base + t * stride + e`,
/// for `t < len` (a power of two) and `e < width`.
#[allow(clippy::too_many_arguments)]
fn fft_rows(
    data: &mut [Complex],
    base: usize,
    len: usize,
    stride: usize,
    width: usize,
    tw: &[Complex],
    rev: &[usize],
    inverse: bool,
) {
    for (i, &j) in rev.iter().enumerate() {
        if i < j {
            for e in 0..width {
                data.swap(base + i * stride + e, base + j * stride + e);
            }
        }
    }
    let mut size = 2;
    while size <= len {
        let half = size / 2;
        let step = len / size;
        for start in (0..len).step_by(size) {
            for k in 0..half {
                let (wr, mut wi) = tw[k * step];
                if inverse {
                    wi = -wi;
                }
                let a0 = base + (start + k) * stride;
                let b0 = base + (start + k + half) * stride;
                for e in 0..width {
                    let (ar, ai) = data[a0 + e];
                    let (br, bi) = data[b0 + e];
                    let (tr, ti) = (br * wr - bi * wi, br * wi + bi * wr);
                    data[a0 + e] = (ar + tr, ai + ti);
                    data[b0 + e] = (ar - tr, ai - ti);
                }
            }
        }
        size <<= 1;
    }
}

/// In-place FFT of a power-of-two length (`inverse` uses `e^{+i...}` and divides by `n`).
pub fn fft(data: &mut [Complex], inverse: bool) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    fft_rows(
        data,
        0,
        n,
        1,
        1,
        &twiddles(n),
        &bit_reverse_index(n),
        inverse,
    );
    if inverse {
        let s = 1.0 / n as f64;
        data.iter_mut().for_each(|v| *v = (v.0 * s, v.1 * s));
    }
}

/// Applies `f` to the first `count` consecutive chunks of length `slab`, in parallel with the
/// `parallel` feature. Each chunk is transformed independently, so results do not depend on
/// the number of threads.
fn for_each_slab(
    data: &mut [Complex],
    slab: usize,
    count: usize,
    f: impl Fn(&mut [Complex]) + Send + Sync,
) {
    let data = &mut data[..slab * count];
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data.par_chunks_mut(slab).for_each(f);
    }
    #[cfg(not(feature = "parallel"))]
    data.chunks_mut(slab).for_each(f);
}

/// FFT along `axis` of a row-major array of shape `n`, only on the lines whose indices on
/// the other axes are below `lim` (no normalisation).
fn fft_axis(data: &mut [Complex], n: [usize; 3], axis: usize, lim: [usize; 3], inverse: bool) {
    let len = n[axis];
    if len <= 1 {
        return;
    }
    let tw = twiddles(len);
    let rev = bit_reverse_index(len);
    match axis {
        // Rows along the last axis, `lim[2]` long; loop over the remaining axis.
        0 => {
            for j in 0..lim[1] {
                fft_rows(data, j * n[2], len, n[1] * n[2], lim[2], &tw, &rev, inverse);
            }
        }
        // One slab per index on axis 0: independent, so they run in parallel.
        1 => for_each_slab(data, n[1] * n[2], lim[0], |slab| {
            fft_rows(slab, 0, len, n[2], lim[2], &tw, &rev, inverse);
        }),
        _ => for_each_slab(data, n[1] * n[2], lim[0], |slab| {
            for j in 0..lim[1] {
                fft_rows(slab, j * n[2], len, 1, 1, &tw, &rev, inverse);
            }
        }),
    }
}

fn check_sizes(n: [usize; 3]) -> Result<()> {
    if n.iter().any(|&s| s > 1 && !s.is_power_of_two()) {
        return invalid(format!("FFT sizes must be powers of two, got {:?}", &n[..]));
    }
    Ok(())
}

/// FFT along every axis of a row-major array of shape `n` (each a power of two or 1).
/// The inverse divides by the number of points.
pub fn fft_nd(data: &mut [Complex], n: [usize; 3], inverse: bool) -> Result<()> {
    fft_nd_pruned(data, n, inverse, n)
}

/// [`fft_nd`] for data that is zero outside the box `[0, nonzero)` (forward), or when only
/// the box `[0, nonzero)` of the result is needed (inverse); values outside the box are then
/// left unspecified.
pub fn fft_nd_pruned(
    data: &mut [Complex],
    n: [usize; 3],
    inverse: bool,
    nonzero: [usize; 3],
) -> Result<()> {
    check_sizes(n)?;
    if inverse {
        // Transform the full grid first, keep only the wanted rows for the later axes.
        fft_axis(data, n, 0, [n[0], n[1], n[2]], true);
        fft_axis(data, n, 1, [nonzero[0], n[1], n[2]], true);
        fft_axis(data, n, 2, [nonzero[0], nonzero[1], n[2]], true);
        let s = 1.0 / (n[0] * n[1] * n[2]) as f64;
        for i in 0..nonzero[0] {
            for j in 0..nonzero[1] {
                let row = (i * n[1] + j) * n[2];
                for v in &mut data[row..row + nonzero[2]] {
                    *v = (v.0 * s, v.1 * s);
                }
            }
        }
    } else {
        // Lines that are entirely zero stay zero: transform the last axis inside the box,
        // then the middle axis for the first `nonzero[0]` slabs, then everything.
        fft_axis(data, n, 2, [nonzero[0], nonzero[1], n[2]], false);
        fft_axis(data, n, 1, [nonzero[0], n[1], n[2]], false);
        fft_axis(data, n, 0, [n[0], n[1], n[2]], false);
    }
    Ok(())
}
