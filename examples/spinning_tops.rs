//! Rigid bodies: the intermediate-axis (Dzhanibekov) flip of a torque-free asymmetric top,
//! and the slow precession of a fast heavy top.
//!
//! Run with `cargo run --release --example spinning_tops`.

use physim::*;

fn main() -> Result<()> {
    // A torque-free top spun almost exactly about its intermediate axis.
    let mut s = RigidSystem::new();
    let mut b = RigidBody::new(1.0, Vec3::new(1.0, 2.0, 3.0));
    b.set_angular_velocity(Vec3::new(0.01, 1.0, 0.02));
    s.add_body(b)?;
    let e0 = s.total_energy()?;
    let traj = s.run(1e-3, 200_000, 1)?;
    let mut flips = Vec::new();
    let mut last = 1.0;
    for f in 0..traj.n_frames() {
        let w = traj.orientation[f]
            .rotate_inverse(traj.angular_velocity[f])
            .y;
        if w.signum() != last {
            flips.push(traj.t[f]);
            last = w.signum();
        }
    }
    let worst = traj
        .total_energy()
        .iter()
        .map(|e| ((e - e0) / e0).abs())
        .fold(0.0, f64::max);
    println!("free asymmetric top, I = (1, 2, 3), spun about axis 2:");
    println!("  intermediate-axis flips at t = {flips:.2?}");
    println!("  worst relative energy error {worst:.1e} (dt = 1e-3, 2e5 steps); L exact");

    // A heavy symmetric top: pivot at the origin, centre of mass at distance l on its axis.
    let (m, g, l, i3, spin) = (1.0, 1.0, 1.0, 0.5, 50.0);
    let mut s = RigidSystem::new();
    let mut top = RigidBody::pivoted(
        m,
        Vec3::new(0.5, 0.5, i3),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, l),
    )?;
    top.orientation = Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.5);
    top.set_angular_velocity(top.to_space(Vec3::new(0.0, 0.0, spin)));
    s.add_body(top)?;
    s.add_force(Box::new(BodyGravity {
        g: Vec3::new(0.0, 0.0, -g),
    }))?;
    let steps = 157_080; // one precession period at the fast-top rate
    let traj = s.run(1e-3, steps, 10)?;
    let mut turned = 0.0;
    for f in 1..traj.n_frames() {
        let a = traj.orientation[f - 1].rotate(Vec3::new(0.0, 0.0, 1.0));
        let b = traj.orientation[f].rotate(Vec3::new(0.0, 0.0, 1.0));
        let mut d = b.y.atan2(b.x) - a.y.atan2(a.x);
        d -= (d / std::f64::consts::TAU).round() * std::f64::consts::TAU;
        turned += d;
    }
    let rate = turned / traj.t[traj.n_frames() - 1];
    let predicted = m * g * l / (i3 * spin);
    println!("heavy symmetric top, spin {spin}, tilt 0.5 rad:");
    println!(
        "  precession rate {rate:.6} vs fast-top mgl/(I3 w3) = {predicted:.6} ({:+.2}%)",
        100.0 * (rate / predicted - 1.0)
    );
    Ok(())
}
