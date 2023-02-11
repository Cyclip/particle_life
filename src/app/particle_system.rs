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

    /// Cell positions
    pub cell_positions: Vec<(usize, usize)>,
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
        let grid = Grid::new(
            constants::GRID_CELL_SIZE as f64,
            (constants::X_SPAN, constants::Y_SPAN),
        );

        let mut behaviour = particle_behaviour::Behaviour::new(num_colours);
        behaviour.randomise();

        let cell_positions = grid.get_cell_positions();

        Particles {
            num_particles: constants::NUM_PARTICLES,
            particles: Array2::zeros((constants::NUM_PARTICLES, 6)),
            grid,
            sim_area,
            graph_area,
            num_colours,
            behaviour,
            cell_positions,
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

        // update every cell
        let cell_positions = self.grid.get_cell_positions();

        self.grid.print();
        
        for cell_pos in cell_positions.iter() {
            println!("Updating cell {:?}", cell_pos);

            self.update_cell(
                cell_pos,
            );

            // print grid
            self.grid.print();
        }
    }

    /// Updates an entire cell, given the particles (and their respective positions)
    /// in its, and its neighbours, cells
    fn update_cell(
        &mut self,
        cell_pos: &(usize, usize),
    ) {
        // get particles in cell
        let particles = self.grid.get_neighbour_particles(cell_pos);

        println!("Cell {:?}: {:?}", cell_pos, particles);

        // O(n^2) algorithm
        // for each particle in cell
        for i in 0..particles.len() {
            // get this particle's pos, vel, colour
            let mut pos = (
                self.particles[[i, 0]],
                self.particles[[i, 1]],
            );

            let mut vel = (
                self.particles[[i, 2]],
                self.particles[[i, 3]],
            );

            // current grid position incase it moves later
            let grid_pos = Grid::sim_pos_to_grid_pos(pos);

            let colour = self.particles[[i, 4]] as usize;

            // for each particle in cell
            for j in 0..particles.len() {
                // if not the same particle
                if i != j {
                    // get this particles pos and colour
                    let other_pos = (
                        self.particles[[j, 0]],
                        self.particles[[j, 1]],
                    );

                    let other_colour = self.particles[[j, 4]] as usize;
                    
                    // get distance between particles
                    let dist = Grid::dist(pos, other_pos);

                    // if particles are close enough
                    let force = if dist < constants::DISTANCE_THRESHOLD {
                        // get force between particles
                        self.behaviour.get_behaviour(colour, other_colour) * constants::FORCE / dist
                    } else {
                        // forced repulsion
                        -1f64 * constants::FORCE / dist
                    };

                    // apply force to velocity
                    vel.0 += force * (other_pos.0 - pos.0) / dist;
                    vel.1 += force * (other_pos.1 - pos.1) / dist;
                }
            }

            // update particle position
            pos.0 += vel.0;
            pos.1 += vel.1;

            // update it in the actual particle array
            self.particles[[i, 0]] = pos.0;
            self.particles[[i, 1]] = pos.1;

            self.particles[[i, 2]] = vel.0;
            self.particles[[i, 3]] = vel.1;

            // wrap around edges
            if pos.0 < 0f64 {
                pos.0 += self.sim_area.0;
            } else if pos.0 > self.sim_area.0 {
                pos.0 -= self.sim_area.0;
            }

            if pos.1 < 0f64 {
                pos.1 += self.sim_area.1;
            } else if pos.1 > self.sim_area.1 {
                pos.1 -= self.sim_area.1;
            }

            // if particle has moved cells
            let new_pos = Grid::sim_pos_to_grid_pos(pos);
            if new_pos != grid_pos {
                self.grid.move_particle(grid_pos, new_pos, i);
            }
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

    /// Get all cell positions
    pub fn get_cell_positions(&self) -> Vec<(usize, usize)> {
        let mut positions = Vec::new();

        for i in 0..self.num_cells.0 {
            for j in 0..self.num_cells.1 {
                positions.push((i, j));
            }
        }

        positions
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
            None => panic!("Particle {} not found in cell: {:?}", particle, pos),
        };

        cell.remove(index);
    }

    /// Moves a particle from one cell to another
    pub fn move_particle(&mut self, old_pos: (usize, usize), new_pos: (usize, usize), particle: usize) {
        println!("Moving particle {} from {:?} to {:?}", particle, old_pos, new_pos);
        self.remove_particle(old_pos, particle);
        self.add_particle(new_pos, particle);
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

    /// Get the particle indexes of all particles in a cell and its neighbours
    pub fn get_neighbour_particles(&self, pos: &(usize, usize)) -> Vec<usize> {
        let mut particles = Vec::new();

        // get cell and neighbours
        let cells = self.get_neighbour_cells_pos(pos);

        // add particles
        for cell in cells {
            particles.append(&mut self.get_cell(cell).clone());
        }

        particles
    }

    /// Get the positions of a cell and its neighbours (including itself)
    fn get_neighbour_cells_pos(&self, pos: &(usize, usize)) -> Vec<(usize, usize)> {
        let mut cells = Vec::new();

        // get cell and neighbours
        let x = pos.0;
        let y = pos.1;

        // prepare edges
        let x_span = self.num_cells.0 - 1;
        let y_span = self.num_cells.1 - 1;

        // add cells
        for i in -1..2 {
            for j in -1..2 {
                // wrap around edges
                let x2 = (
                    x as i32 + i + x_span as i32
                ) as usize % x_span;

                let y2 = (
                    y as i32 + j + y_span as i32
                ) as usize % y_span;

                cells.push((x2, y2));
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

    /// Returns the distance between two simulation positions
    pub fn dist(pos1: (f64, f64), pos2: (f64, f64)) -> f64 {
        let x = (pos1.0 as f64 - pos2.0 as f64).powi(2);
        let y = (pos1.1 as f64 - pos2.1 as f64).powi(2);

        (x + y).sqrt()
    }

    /// Prints the grid in a grid fashion
    pub fn print(&self) {
        // 6 cells per line
        let mut line = String::new();

        for i in 0..self.num_cells.0 {
            for j in 0..self.num_cells.1 {
                let cell = self.get_cell((i, j));
                // list of all particles in cell
                let mut cell_str = String::new();

                for particle in cell {
                    cell_str.push_str(&format!("{}, ", particle));
                }

                line.push_str(&format!("| {:<20} ", cell_str));
            }

            println!("{}", line);
            line.clear();
        }

        println!();
    }
}