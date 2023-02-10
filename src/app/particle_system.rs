//! Handles particle system

use crate::constants;
use crate::app::particle_behaviour;

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

    /// Behaviour of particles
    pub behaviour: particle_behaviour::Behaviour,
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
        let mut behaviour = particle_behaviour::Behaviour::new(num_colours);
        behaviour.randomise();

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
            behaviour,
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

            // get grid cell index
            let grid_pos = Grid::sim_pos_to_grid_pos((x, y));

            self.particles[[i, 0]] = x;
            self.particles[[i, 1]] = y;

            // add particle to grid cell
            self.grid.add_particle(grid_pos, i);
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
        // apply resistance to all particles (velocity)
        for i in 0..self.num_particles {
            self.particles[[i, 2]] *= constants::RESISTANCE;
            self.particles[[i, 3]] *= constants::RESISTANCE;
        }

        // to store particles that have moved cells
        // (index, (old cell), (new cell))
        let mut moved_particles: Vec<(usize, (usize, usize), (usize, usize))> = Vec::new();

        // iterate through all the grids cells
        for (pos, _) in self.grid.cells.iter_mut() {
            // get all particles in current and neighbour cells
            let particles = self.grid.get_neighbour_cells(pos);

            for i in 0..particles.len() {
                // get particle from index
                let particle = particles[i];

                // get particle type
                let particle_type = self.particles[[particle, 4]] as usize;

                // get particle position
                let particle_pos = (self.particles[[particle, 0]], self.particles[[particle, 1]]);

                // get particle velocity
                let particle_vel = (self.particles[[particle, 2]], self.particles[[particle, 3]]);

                // iterate through all other particles in cell
                for j in 0..particles.len() {
                    if i == j {
                        continue;
                    } else {
                        // get particle from index
                        let other_particle = particles[j];

                        // get particle type
                        let other_particle_type = self.particles[[other_particle, 4]] as usize;

                        // get particle position
                        let other_particle_pos = (
                            self.particles[[other_particle, 0]],
                            self.particles[[other_particle, 1]],
                        );

                        // get resulting behaviour
                        let resulting_behaviour = self.behaviour.get_behaviour(
                            particle_type,
                            other_particle_type,
                        );

                        // apply behaviour
                        // add resulting_behaviour * distance_dif * FORCE to velocity
                        let distance_dif = (
                            other_particle_pos.0 - particle_pos.0,
                            other_particle_pos.1 - particle_pos.1,
                        );

                        let distance = (distance_dif.0.powi(2) + distance_dif.1.powi(2)).sqrt();

                        let force = resulting_behaviour * distance_dif.0 * constants::FORCE / distance;
                        
                        self.particles[[particle, 2]] += force;
                        self.particles[[particle, 3]] += force;
                    }
                }

                // change position
                self.particles[[particle, 0]] += particle_vel.0;
                self.particles[[particle, 1]] += particle_vel.1;

                // wrap around edges
                if self.particles[[particle, 0]] < 0.0 {
                    self.particles[[particle, 0]] += self.sim_area.0;
                } else if self.particles[[particle, 0]] > self.sim_area.0 {
                    self.particles[[particle, 0]] -= self.sim_area.0;
                }

                if self.particles[[particle, 1]] < 0.0 {
                    self.particles[[particle, 1]] += self.sim_area.1;
                } else if self.particles[[particle, 1]] > self.sim_area.1 {
                    self.particles[[particle, 1]] -= self.sim_area.1;
                }

                // get new grid position
                let new_grid_pos = Grid::sim_pos_to_grid_pos(particle_pos);

                // check if particle has moved cells
                if new_grid_pos != *pos {
                    moved_particles.push((particle, *pos, new_grid_pos));
                }
            }
        }

        // move particles to new grid cells
        for (particle, old_pos, new_pos) in moved_particles {
            self.grid.remove_particle(old_pos, particle);
            self.grid.add_particle(new_pos, particle);

            println!("Moved particle {} from {:?} to {:?}", particle, old_pos, new_pos);
        }
    }
}

impl Grid {
    pub fn new(cell_size: f64, area: (f64, f64)) -> Grid {
        let num_cells = (
            (area.0 / cell_size).ceil() as usize,
            (area.1 / cell_size).ceil() as usize,
        );

        let mut cells = HashMap::new();

        let start = num_cells.0 + 1;
        let end = num_cells.1 + 1;

        for i in 0..start {
            for j in 0..end {
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

    /// Adds a particle to a grid cell
    pub fn add_particle(&mut self, pos: (usize, usize), particle: usize) {
        // if not already in cell
        if !self.get_cell(pos).contains(&particle) {
            self.get_cell_mut(pos).push(particle);
        }
    }

    /// Removes a particle from a grid cell
    pub fn remove_particle(&mut self, pos: (usize, usize), particle: usize) {
        let cell = self.get_cell_mut(pos);

        let index = match cell.iter().position(|&x| x == particle) {
            Some(index) => index,
            None => panic!("Particle not found in cell: {:?}", pos),
        };

        cell.remove(index);
    }

    /// Returns a reference to a grid cell
    pub fn get_cell(&self, pos: (usize, usize)) -> &Vec<usize> {
        match self.cells.get(&pos) {
            Some(cell) => cell,
            None => panic!("Cell not found: {:?}", pos),
        }
    }

    /// Returns a mutable reference to a grid cell
    pub fn get_cell_mut(&mut self, pos: (usize, usize)) -> &mut Vec<usize> {
        match self.cells.get_mut(&pos) {
            Some(cell) => cell,
            None => panic!("Cell not found: {:?}", pos),
        }
    }

    /// Get the positions and particles within a cell and its neighbours (including itself)
    /// Returns a vector of format:
    /// [(pos, particles), (pos, particles), ...]
    pub fn get_neighbour_cells(&self, pos: &(usize, usize)) -> Vec<((usize, usize), &Vec<usize>)> {
        let mut cells = Vec::new();

        let start = (pos.0 - 1, pos.1 - 1);
        let end = (pos.0 + 1, pos.1 + 1);

        for i in start.0..=end.0 {
            for j in start.1..=end.1 {
                let cell_pos = (i, j);

                // if cell is within grid
                if cell_pos.0 <= self.num_cells.0 && cell_pos.1 <= self.num_cells.1 {
                    cells.push((cell_pos, self.get_cell(cell_pos)));
                }
            }
        }

        cells
    }

    /// Returns grid cell index of a given simulation position
    pub fn sim_pos_to_grid_pos(pos: (f64, f64)) -> (usize, usize) {
        let x = (pos.0 / constants::GRID_CELL_SIZE as f64).floor() as usize;
        let y = (pos.1 / constants::GRID_CELL_SIZE as f64).floor() as usize;

        (x, y)
    }

    /// Returns simulation position of a given grid cell index (approximate)
    pub fn grid_pos_to_sim_pos(pos: (usize, usize)) -> (f64, f64) {
        let x = pos.0 as f64 * constants::GRID_CELL_SIZE as f64;
        let y = pos.1 as f64 * constants::GRID_CELL_SIZE as f64;

        (x, y)
    }
}