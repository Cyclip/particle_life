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

const BLACK: [f32; 4] = [0.0, 0.0, 0.0, 1.0];

/// Controls the graphics aspect of the program
/// Handles updates and rendering
pub struct App {
    pub gl: GlGraphics, // OpenGL drawing backend
}

impl App {
    /// Creates a new App
    pub fn new(opengl: OpenGL) -> App {
        App {
            gl: GlGraphics::new(opengl),
        }
    }

    /// Renders a frame
    pub fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

    }

    /// Updates the state of the App
    pub fn update(&mut self, args: &UpdateArgs) {
    }
}