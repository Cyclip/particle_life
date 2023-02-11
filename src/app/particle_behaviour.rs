//! Manages the behaviour of particles based on their colour

use crate::constants;

use ndarray::Array2;

/// Represents all the behaviours of particles
/// Stores a 2D matrix of behaviours of (nxn) where n is the number of colours
pub struct Behaviour {
    pub behaviours: Array2<f64>,
}

impl Behaviour {
    pub fn new(num_colours: usize) -> Behaviour {
        // create a zero matrix of size (num_colours x num_colours)
        let behaviours = Array2::zeros((num_colours, num_colours));
        
        Behaviour { behaviours }
    }

    /// Randomise all behaviours from -1 to 1
    pub fn randomise(&mut self) {
        for i in 0..self.behaviours.shape()[0] {
            for j in 0..self.behaviours.shape()[1] {
                self.behaviours[[i, j]] = rand::random::<f64>() * 2.0 - 1.0;
            }
        }

        // pretty print
        println!("{:#?}", self.behaviours);
    }

    /// Get the behaviour of a particle with colour `colour` towards a particle with colour `target_colour`
    pub fn get_behaviour(&self, colour: usize, target_colour: usize) -> f64 {
        self.behaviours[[colour, target_colour]]
    }

    /// Set the behaviour of a particle with colour `colour` towards a particle with colour `target_colour`
    pub fn set_behaviour(&mut self, colour: usize, target_colour: usize, behaviour: f64) {
        self.behaviours[[colour, target_colour]] = behaviour;
    }

    /// Get all behaviours of a particle with colour `colour`
    pub fn get_behaviours(&self, colour: usize) -> Vec<f64> {
        let mut behaviours = Vec::new();

        for i in 0..self.behaviours.shape()[1] {
            behaviours.push(self.behaviours[[colour, i]]);
        }

        behaviours
    }

    /// Get the new velocity of a particle with colour `colour` towards a particle with colour `target_colour`
    pub fn get_new_velocity(
        &self,
        behaviour: f64,
        velocity: (f64, f64),
    ) -> (f64, f64) {

        let new_velocity = (
            velocity.0 + behaviour,
            velocity.1 + behaviour,
        );

        new_velocity
    }
}