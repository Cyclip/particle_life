//! Stores constants used in the program

// ================== Simulation ==================
/// Number of particles in the simulation
pub const NUM_PARTICLES: usize = 500;
/// x-span of the simulation
pub const X_SPAN: f64 = 100.0;
/// y-span of the simulation
pub const Y_SPAN: f64 = 100.0;
/// Resistance multiplier (friction)
pub const RESISTANCE: f64 = 0.9;
/// Force multiplier
pub const FORCE: f64 = 0.0001f64;
/// Number of colours
pub const NUM_COLOURS: usize = 10;
/// Spatial partitioning grid size (must be a divisor of X_SPAN and Y_SPAN)
pub const GRID_CELL_SIZE: usize = 10;
/// All available colours
pub const COLOURS: [[f32; 4]; 12] = [
    [1.0, 0.0, 0.0, 1.0], // red
    [0.0, 1.0, 0.0, 1.0], // green
    [0.0, 0.0, 1.0, 1.0], // blue
    [1.0, 1.0, 0.0, 1.0], // yellow
    [1.0, 0.0, 1.0, 1.0], // magenta
    [0.0, 1.0, 1.0, 1.0], // cyan
    [1.0, 0.5, 0.0, 1.0], // orange
    [0.5, 0.0, 1.0, 1.0], // purple
    [0.0, 1.0, 0.5, 1.0], // teal
    [0.5, 1.0, 0.0, 1.0], // lime
    [1.0, 0.0, 0.5, 1.0], // pink
    [0.0, 0.5, 1.0, 1.0], // sky blue
];
/// Particle size
pub const PARTICLE_SIZE: f64 = 2.0;
/// Maximum distance between particles before forced repulsion
pub const DISTANCE_THRESHOLD: f64 = 1.0;
// ================== Simulation ==================


// ================== Window ==================
/// Width of the window
pub const WIDTH: u32 = 1200;
/// Height of the window
pub const HEIGHT: u32 = 700;
// ================== Window ==================