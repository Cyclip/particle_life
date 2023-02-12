#[allow(unused_imports)]
use glutin_window::GlutinWindow as Window;
#[allow(unused_imports)]
use opengl_graphics::{GlGraphics, OpenGL};
#[allow(unused_imports)]
use piston::event_loop::{EventSettings, Events};
#[allow(unused_imports)]
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent, MouseCursorEvent};
#[allow(unused_imports)]
use piston::window::WindowSettings;

use crate::app::particle_system::Grid;
use crate::constants;

mod particle_system;
pub mod particle_behaviour;

const BLACK: [f32; 4] = [0.0, 0.0, 0.0, 1.0];

/// Controls the graphics aspect of the program
/// Handles updates and rendering
pub struct App {
    pub gl: GlGraphics, // OpenGL drawing backend
    particles: particle_system::Particles,

    last_mouse_pos: (f64, f64),
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

        App {
            gl: GlGraphics::new(opengl),
            particles,
            last_mouse_pos: (0f64, 0f64),
        }
    }

    /// Mouse cursor event
    pub fn mouse_cursor(&mut self, args: &[f64; 2]) {
        // update
        self.last_mouse_pos = (args[0], args[1]);
    }

    /// Renders a frame
    pub fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        // draw all particles
        self.gl.draw(args.viewport(), |c, gl| {
            clear(BLACK, gl);

            let window_width = args.window_size[0];
            let window_height = args.window_size[1];

            // draw currently selected grid
            // screen -> sim -> grid -> sim -> screen
            
            // screen -> sim
            let mouse_sim_pos = App::screen_pos_to_sim_pos(
                self.last_mouse_pos, 
                (window_width, window_height)
            );

            // sim -> grid
            let mouse_grid_pos = Grid::sim_pos_to_grid_pos(mouse_sim_pos);

            // grid -> sim
            let mouse_sim_pos = Grid::grid_pos_to_sim_pos(mouse_grid_pos);

            // sim -> screen
            let mouse_screen_pos = App::sim_pos_to_screen_pos(mouse_sim_pos, (window_width, window_height));

            let size = App::sim_pos_to_screen_pos(
                (
                    constants::GRID_CELL_SIZE as f64,
                    constants::GRID_CELL_SIZE as f64,
                ),
                (window_width, window_height),
            );

            rectangle(
                [1.0, 1.0, 1.0, 0.01],
                [
                    mouse_screen_pos.0 - constants::GRID_CELL_SIZE as f64 / 2f64,
                    mouse_screen_pos.1 - constants::GRID_CELL_SIZE as f64 / 2f64,
                    size.0,
                    size.1,
                ],
                c.transform,
                gl,
            );


            // draw all particles
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
        // update particles
        self.particles.update();
    }

    /// Converts a position in the simulation to a position on the screen
    fn sim_pos_to_screen_pos(pos: (f64, f64), screen_size: (f64, f64)) -> (f64, f64) {
        let x = (pos.0 / constants::X_SPAN) * screen_size.0;
        let y = (pos.1 / constants::Y_SPAN) * screen_size.1;

        (x + 2f64, y + 2f64)
    }

    /// Converts a position on the screen to a position in the simulation
    fn screen_pos_to_sim_pos(pos: (f64, f64), screen_size: (f64, f64)) -> (f64, f64) {
        let x = (pos.0 / screen_size.0) * constants::X_SPAN;
        let y = (pos.1 / screen_size.1) * constants::Y_SPAN;

        (x, y)
    }
}