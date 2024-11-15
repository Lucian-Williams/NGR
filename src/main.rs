#![allow(dead_code)]

extern crate glutin_window;
extern crate graphics;
extern crate opengl_graphics;
extern crate piston;

use crate::physics::*;

use libm::pow;
use std::env;
use std::collections::LinkedList;
use array_init::array_init;
use glutin_window::GlutinWindow as Window;
use opengl_graphics::{GlGraphics, OpenGL};
use piston::event_loop::{EventSettings, Events};
use piston::input::{RenderArgs, RenderEvent, UpdateArgs, UpdateEvent};
use piston::window::WindowSettings;
use graphics::math::Vec2d;

mod vector2d;
mod physics;

const NUM_PARTICLES: usize = 1024;
const NUM_POINTS: usize = 16;
const SPACE_SCALE: f64 = 0.01;
const TIME_STEP: f64 = 0.0000001;
const STEPS_PER_FRAME: usize = 1;
const STEPS_PER_LINE: usize = 64;
const R0: f64 = 2.0 * GM / C2;

pub struct App {
    gl: GlGraphics,
    particles: [(Particle, LinkedList<Vec2d>); NUM_PARTICLES],
    step_count: usize,
    impact_count: u32 // Number of phtons that have hit the event horizon
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        const WHITE: [f32; 4] = [1.0; 4];
        const BLACK: [f32; 4] = [0.0, 0.0, 0.0, 1.0];

        let (x, y) = (args.window_size[0] / 2.0, args.window_size[1] / 2.0);

        self.gl.draw(args.viewport(), |c, gl| {
            clear(BLACK, gl);

            let transform = c
                .transform
                .trans(x, y);

            for particle in &self.particles {
                let mut yellow: [f32; 4] = [1.0, 1.0, 0.0, 1.0];
                if let Some(front) = particle.1.front() {
                    let mut prev = front;
                    for (i, point) in particle.1.iter().enumerate() {
                        if i != 0 {
                            line_from_to(
                                yellow,
                                1.0,
                                [prev[0], prev[1]],
                                [point[0], point[1]],
                                transform,
                                gl
                            );

                            prev = point;
                            yellow[0] *= 0.7;
                            yellow[1] *= 0.7;
                        }
                    }
                }
            }

            ellipse_from_to(
                    WHITE,
                    [-SPACE_SCALE * 2.0 * GM / C2; 2],
                    [SPACE_SCALE * 2.0 * GM / C2; 2],
                    transform, gl
            );
        });
    }

    fn update(&mut self, _args: &UpdateArgs) {
        let mut inside_count: u32 = 0;

        for particle in &mut self.particles {
            particle.0.update(TIME_STEP);
            particle.0.correct_velocity();

            if let Some(front) = particle.1.front_mut() {
                *front = [
                    SPACE_SCALE * particle.0.r.x,
                    SPACE_SCALE * particle.0.r.y
                ];
            }

            if self.step_count % STEPS_PER_LINE == 0 {
                particle.1.push_front([
                    SPACE_SCALE * particle.0.r.x,
                    SPACE_SCALE * particle.0.r.y,
                ]);
                if particle.1.len() > NUM_POINTS {
                    particle.1.pop_back();
                }
            }

            if particle.0.r.sqr_magnitude() <= R0 * R0 {
                inside_count += 1;
            }
        }

        if inside_count < self.impact_count {
            println!("\noops");
        }

        else {
            self.impact_count = inside_count;
            print!("\x08\x08\x08\x08{}", inside_count);
        }

        self.step_count += 1;
    }
}

fn newtons_method<F>(f: F, mut x: f64, mut n: i32) -> f64 where
    F: Fn(f64) -> f64 {

    while n > 0 {
        x = x - f(x);
        n -= 1;
    }

    x
}

fn first_guess(x: f64) -> f64 {
    pow(x - 0.5, 5.0) + x / 4.0 + 3.0 / 8.0
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut span: f64 = 20.0;
    let mut proximity: f64 = 20.0;
    if args.len() > 1 {
        if let Ok(i) = args[1].parse::<f64>() {
            span = i;
        }
    }
    if args.len() > 2 {
        if let Ok(i) = args[2].parse::<f64>() {
            proximity = i;
        }
    }

    fn u(x: f64) -> f64 {
        pow(x, 5.0) - pow(x - 1.0, 5.0)
    }
    let mut c: f64 = 3.0 / span;
    let mut step: f64 = 0.0625;

    for _i in 0..60 {
        let mut particle: Particle = Particle::new_coord_launch(proximity * R0, c * -span * R0, -1.0, 0.0);
        while particle.r.sqr_magnitude() > 2.0 * R0 * R0 && particle.v.x * particle.r.x + particle.v.y * particle.r.y < 0.0 {
            particle.update(TIME_STEP);
            particle.correct_velocity();
        }
        if particle.r.sqr_magnitude() <= 2.0 * R0 * R0 {
            c += step;
        }
        else {
            c -= step;
        }
        step *= 0.5;
    }

    let b: f64 = newtons_method(|x| {
            (pow(x, 5.0) - c * (u(x))) * (u(x)) /
            (5.0 * pow(x, 4.0) * pow(x - 1.0, 4.0))
            }, first_guess(c), 10);
    let a: f64 = c / pow(b, 5.0);
    let opengl = OpenGL::V3_2;
    
    let mut window: Window = WindowSettings::new("Light passing by the Schwarzschild black hole", [1920, 1080])
        .graphics_api(opengl)
        .exit_on_esc(true)
        //.fullscreen(true)
        .build()
        .unwrap();

    let particle_array: [(Particle, LinkedList<Vec2d>); NUM_PARTICLES] = array_init(|i| {
        let particle: Particle = Particle::new_coord_launch(
            proximity * R0,
            (a * pow(i as f64 / (NUM_PARTICLES - 1 ) as f64 - b, 5.0) + c) * -span * R0,
            //i as f64 / (NUM_PARTICLES - 1 ) as f64 * -4.0 * R0,
            -1.0,
            0.0
        );
        let mut init_points: LinkedList<Vec2d> = LinkedList::new();
        for _i in 0..2 {
            init_points.push_front([
                SPACE_SCALE * particle.r.x,
                SPACE_SCALE * particle.r.y
            ]);
        }
        (particle, init_points)
    });

    println!("Impact Count:");

    let mut app = App {
        gl: GlGraphics::new(opengl),
        particles: particle_array,
        step_count: 0,
        impact_count: 0
    };

    let mut events = Events::new(EventSettings::new());
    while let Some(e) = events.next(&mut window) {
        if let Some(args) = e.render_args() {
            app.render(&args);
        }

        if let Some(args) = e.update_args() {
            app.update(&args);
        }
    }
}
