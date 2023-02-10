extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;
extern crate ndarray;

mod app;
pub mod constants;

use glutin_window::GlutinWindow as Window;
#[allow(unused_imports)]
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
#[allow(unused_imports)]
use piston::input::{RenderEvent, UpdateEvent, RenderArgs, UpdateArgs, MouseCursorEvent};
use piston::window::WindowSettings;

use app::App;

fn main() {
    let opengl = OpenGL::V3_2;

    // Create a Glutin window.
    let mut window: Window = WindowSettings::new(
        "Particle life", 
        [constants::WIDTH, constants::HEIGHT]
    )
        .graphics_api(opengl)
        .exit_on_esc(true)
        .build()
        .unwrap();

    // Create a new game and run it.
    let mut app = App::new(opengl);

    let mut events = Events::new(EventSettings::new());

    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args);
        }

        if let Some(args) = e.mouse_cursor_args() {
            app.mouse_cursor(&args);
        }
    }
}