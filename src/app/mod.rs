#[allow(unused_imports)]
use glutin_window::GlutinWindow as Window;
#[allow(unused_imports)]
use opengl_graphics::{GlGraphics, OpenGL};
#[allow(unused_imports)]
use piston::event_loop::{EventSettings, Events};
#[allow(unused_imports)]
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
#[allow(unused_imports)]
use piston::window::WindowSettings;

use crate::constants;

mod particle_system;

const BLACK: [f32; 4] = [0.0, 0.0, 0.0, 1.0];

/// Controls the graphics aspect of the program
/// Handles updates and rendering
pub struct App {
    pub gl: GlGraphics, // OpenGL drawing backend
    particles: particle_system::Particles,
}

impl App {
    /// Creates a new App
    pub fn new(opengl: OpenGL) -> App {
        let mut particles: particle_system::Particles = particle_system::Particles::new(
            (constants::X_SPAN, constants::Y_SPAN),
            (constants::WIDTH as f64, constants::HEIGHT as f64),
            constants::NUM_COLOURS,
        );
        particles.init();

        // print pos with pos and colour
        for i in 0..particles.num_particles {
            println!(
                "pos: ({}, {}), colour: {}",
                particles.particles[[i, 0]],
                particles.particles[[i, 1]],
                particles.particles[[i, 4]]
            );
        }

        App {
            gl: GlGraphics::new(opengl),
            particles,
        }
    }

    /// Renders a frame
    pub fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        // draw all particles
        self.gl.draw(args.viewport(), |c, gl| {
            clear(BLACK, gl);

            let window_width = args.window_size[0];
            let window_height = args.window_size[1];

            for i in 0..self.particles.num_particles {
                // convert simulation position to graphical position
                let pos = App::sim_pos_to_screen_pos((
                    self.particles.particles[[i, 0]],
                    self.particles.particles[[i, 1]],
                ), 
                (window_width, window_height
            ));

                let colour = constants::COLOURS[self.particles.particles[[i, 4]] as usize];

                // size depending on screen size
                let size = constants::PARTICLE_SIZE * (window_width / constants::WIDTH as f64) * 2f64;

                // draw a circle
                ellipse(
                    colour,
                    [
                        pos.0 - constants::PARTICLE_SIZE,
                        pos.1 - constants::PARTICLE_SIZE,
                        size,
                        size,
                    ],
                    c.transform,
                    gl,
                );
            }
        });
    }

    /// Updates the state of the App
    pub fn update(&mut self, args: &UpdateArgs) {
    }

    /// Converts a position in the simulation to a position on the screen
    fn sim_pos_to_screen_pos(pos: (f64, f64), screen_size: (f64, f64)) -> (f64, f64) {
        let x = (pos.0 / constants::X_SPAN) * screen_size.0;
        let y = (pos.1 / constants::Y_SPAN) * screen_size.1;

        

        (x, y)
    }
}