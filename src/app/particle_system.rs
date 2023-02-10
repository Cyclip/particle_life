//! Handles particle system

use crate::constants;

use ndarray::Array2;
use std::collections::HashMap;

/// Stores all particles within an array.
/// All particles contain the following properties:
/// 0, 1:   Position          (x, y)
/// 2, 3:   Velocity          (vx, vy)
/// 4:      Type (colour)     (int)
/// 
/// These are stored within a single array (for performance).
/// Here is how the array is structured:
/// [
///     x1, y1, vx1, vy1, type1
///     x2, y2, vx2, vy2, type2
///     ...
/// ]
/// 
/// Particles' index in the array is the same as their index in the hashmap.
pub struct Particles {
    /// Number of particles
    pub num_particles: usize,

    /// Array of particles
    pub particles: Array2<f64>,

    /// Hashmap of grid cells
    pub grid: Grid,

    /// Simulation area
    pub sim_area: (f64, f64),

    /// Graphical area
    pub graph_area: (f64, f64),

    /// Number of colours
    pub num_colours: usize,
}

/// Hashmap of grid cells
/// Each cell contains a vector of particle indices
pub struct Grid {
    /// Hashmap of grid cells (2d)
    pub cells: HashMap<(usize, usize), Vec<usize>>,

    /// Size of each grid cell
    pub cell_size: f64,

    /// Simulation area
    pub area: (f64, f64),

    /// Number of grid cells in x and y direction
    pub num_cells: (usize, usize),
}

impl Particles {
    pub fn new(sim_area: (f64, f64), graph_area: (f64, f64), num_colours: usize) -> Particles {
        Particles {
            num_particles: constants::NUM_PARTICLES,
            particles: Array2::zeros((constants::NUM_PARTICLES, 6)),
            grid: Grid::new(
                constants::GRID_CELL_SIZE as f64,
                (constants::X_SPAN, constants::Y_SPAN),
            ),
            sim_area,
            graph_area,
            num_colours,
        }
    }

    /// Initialises all particles
    pub fn init(&mut self) {
        self.random_positions();
        self.random_colours();
    }

    /// Sets all particles to random positions within the simulation area
    pub fn random_positions(&mut self) {
        for i in 0..self.num_particles {
            let x = rand::random::<f64>() * self.sim_area.0;
            let y = rand::random::<f64>() * self.sim_area.1;

            self.particles[[i, 0]] = x;
            self.particles[[i, 1]] = y;
        }
    }

    /// Sets all particles to random colours
    pub fn random_colours(&mut self) {
        for i in 0..self.num_particles {
            let colour = rand::random::<usize>() % self.num_colours;

            self.particles[[i, 4]] = colour as f64;
        }
    }

    /// Update all grid cells
    pub fn update(&mut self) {

    }
}

impl Grid {
    pub fn new(cell_size: f64, area: (f64, f64)) -> Grid {
        let num_cells = (
            (area.0 / cell_size).ceil() as usize,
            (area.1 / cell_size).ceil() as usize,
        );

        let mut cells = HashMap::new();

        for i in 0..num_cells.0 {
            for j in 0..num_cells.1 {
                cells.insert((i, j), Vec::new());
            }
        }

        Grid {
            cells,
            cell_size,
            area,
            num_cells,
        }
    }

    /// Returns a reference to a grid cell
    pub fn get_cell(&self, pos: (usize, usize)) -> &Vec<usize> {
        self.cells.get(&pos).unwrap()
    }

    /// Returns a mutable reference to a grid cell
    pub fn get_cell_mut(&mut self, pos: (usize, usize)) -> &mut Vec<usize> {
        self.cells.get_mut(&pos).unwrap()
    }

    /// Get all particles in cell pos
    pub fn get_particles(&self, pos: (usize, usize)) -> Vec<usize> {
        self.get_cell(pos).clone()
    }

    /// Get all particles in cell pos and neighbouring cells
    pub fn get_particles_neighbours(&self, pos: (usize, usize)) -> Vec<usize> {
        let mut particles = Vec::new();

        for i in pos.0 - 1..pos.0 + 2 {
            for j in pos.1 - 1..pos.1 + 2 {
                if i < self.num_cells.0 && j >= 0 && j < self.num_cells.1 {
                    particles.append(&mut self.get_particles((i, j)).clone());
                }
            }
        }

        particles
    }

    /// Returns grid cell index of a given simulation position
    pub fn sim_pos_to_grid_pos(&self, pos: (f64, f64)) -> (usize, usize) {
        let x = (pos.0 / self.cell_size).floor() as usize;
        let y = (pos.1 / self.cell_size).floor() as usize;

        (x, y)
    }

    /// Returns simulation position of a given grid cell index (approximate)
    pub fn grid_pos_to_sim_pos(&self, pos: (usize, usize)) -> (f64, f64) {
        let x = pos.0 as f64 * self.cell_size;
        let y = pos.1 as f64 * self.cell_size;

        (x, y)
    }
}