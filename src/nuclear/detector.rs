//! Detector faces for the per-event loop: the same geometry as `physim.nuclear.detectors.Geometry`
//! (Python), which passes in its centre and axes so both sides segment faces identically.

use crate::vec3::Vec3;

/// Face shape and segmentation.
#[derive(Clone, Debug)]
pub enum Shape {
    /// `width` × `height` (mm), `strips_x` along u and `strips_y` along v.
    Rectangle {
        width: f64,
        height: f64,
        strips_x: u32,
        strips_y: u32,
    },
    /// Annulus (a disc when `inner = 0`), `rings` of equal radial width and `sectors` of equal
    /// angle counted from u towards v.
    Annular {
        inner: f64,
        outer: f64,
        rings: u32,
        sectors: u32,
    },
}

/// One detector: a flat face with its response parameters.
#[derive(Clone, Debug)]
pub struct Face {
    pub centre: Vec3,
    /// Unit normal pointing from the face towards the target.
    pub n: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub shape: Shape,
    /// Index of the detector material in the generator's material list.
    pub material: usize,
    /// Dead layer and active thickness, mg/cm².
    pub dead_layer: f64,
    pub thickness: f64,
    /// Resolution (FWHM) and threshold, MeV.
    pub fwhm: f64,
    pub threshold: f64,
}

/// Where a track crosses a face.
#[derive(Clone, Copy, Debug)]
pub struct Hit {
    pub distance: f64,
    pub segment: (u32, u32),
    /// cos of the angle between the track and the face normal.
    pub cos_incidence: f64,
}

impl Face {
    /// Crossing of the line from `source` along unit vector `d` with the front of the face.
    pub fn hit(&self, source: Vec3, d: Vec3) -> Option<Hit> {
        let dn = d.dot(self.n);
        if dn >= 0.0 {
            return None; // parallel, or arriving from behind
        }
        let t = (self.centre - source).dot(self.n) / dn;
        if t <= 0.0 {
            return None;
        }
        let rel = source + d * t - self.centre;
        let (su, sv) = (rel.dot(self.u), rel.dot(self.v));
        let segment = match self.shape {
            Shape::Rectangle {
                width,
                height,
                strips_x,
                strips_y,
            } => {
                let (a, b) = (width / 2.0, height / 2.0);
                if su.abs() > a || sv.abs() > b {
                    return None;
                }
                (
                    bin((su + a) / width, strips_x),
                    bin((sv + b) / height, strips_y),
                )
            }
            Shape::Annular {
                inner,
                outer,
                rings,
                sectors,
            } => {
                let r = su.hypot(sv);
                if r < inner || r > outer {
                    return None;
                }
                let ang = sv.atan2(su).rem_euclid(std::f64::consts::TAU);
                (
                    bin((r - inner) / (outer - inner), rings),
                    bin(ang / std::f64::consts::TAU, sectors),
                )
            }
        };
        Some(Hit {
            distance: t,
            segment,
            cos_incidence: -dn,
        })
    }
}

/// Which of `n` equal bins the fraction `f` in [0, 1] falls in (clamped, as in the Python code).
fn bin(f: f64, n: u32) -> u32 {
    ((f * n as f64) as i64).clamp(0, n as i64 - 1) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit(v: Vec3) -> Vec3 {
        v / v.norm()
    }

    #[test]
    fn strips_and_rings() {
        let rect = Face {
            centre: Vec3::new(0.0, 0.0, 100.0),
            n: Vec3::new(0.0, 0.0, -1.0),
            u: Vec3::new(-1.0, 0.0, 0.0),
            v: Vec3::new(0.0, 1.0, 0.0),
            shape: Shape::Rectangle {
                width: 40.0,
                height: 20.0,
                strips_x: 4,
                strips_y: 2,
            },
            material: 0,
            dead_layer: 0.0,
            thickness: 1.0,
            fwhm: 0.0,
            threshold: 0.0,
        };
        let h = rect
            .hit(Vec3::ZERO, unit(Vec3::new(0.15, 0.05, 1.0)))
            .unwrap();
        // u = −x: x = 15 mm → su = −15 → strip 0; y = 5 → sv = 5 → strip 1.
        assert_eq!(h.segment, (0, 1));
        assert!(
            (h.distance - (0.15f64.powi(2) + 0.05f64.powi(2) + 1.0).sqrt() * 100.0).abs() < 1e-9
        );
        assert!(rect.hit(Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0)).is_none());
        assert!(rect
            .hit(Vec3::ZERO, unit(Vec3::new(0.3, 0.0, 1.0)))
            .is_none());
        let mut ring = rect.clone();
        ring.u = Vec3::new(1.0, 0.0, 0.0);
        ring.shape = Shape::Annular {
            inner: 10.0,
            outer: 30.0,
            rings: 4,
            sectors: 8,
        };
        // Radius 27 mm at 100 degrees: ring 3, sector 2 (90-135 degrees).
        let a = 100f64.to_radians();
        let d = unit(Vec3::new(27.0 * a.cos(), 27.0 * a.sin(), 100.0));
        assert_eq!(ring.hit(Vec3::ZERO, d).unwrap().segment, (3, 2));
        assert!(ring.hit(Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)).is_none());
    }
}
