//! The nuclear experiment planner's per-event engine: two-body kinematics, energy loss from
//! tables, detector faces and the Monte Carlo event generator. The physics models themselves
//! (cross sections, stopping powers, straggling) live in the Python package
//! `physim.nuclear`, which tabulates them and hands the tables to [`events::Generator`].

pub mod detector;
pub mod events;
pub mod kinematics;
pub mod stopping;

pub use detector::{Face, Shape};
pub use events::{Channel, Generator, Layer, Record};
pub use kinematics::TwoBody;
pub use stopping::Table;
