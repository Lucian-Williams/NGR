use libm::sqrt;
use libm::pow;
use crate::vector2d::Vector2D;

// Currently, support for massive particle behavior is not available for
// anisotropic lightspeed, so all massive particle initialization is disabled.

pub const G: f64 = 0.000000000066743;                   // Gravitational constant
pub const M: f64 = 1989100000000000000000000000000.0;   // Sun's mass
pub const GM: f64  = G * M;                             // GM
pub const C: f64 = 299792458.0;                         // c
pub const C2: f64 = C * C;                              // c^2

// Schwarschild time dilation squared ("alpha sub g",
// where 'g' is for "gravity")
fn sqr_alpha_g(r: &Vector2D) -> f64 {
    if r.magnitude() == 0.0 {return 0.0;}

    let val = 1.0 - 2.0 * GM / (r.magnitude() * C2);
    if val > 0.0 {val} else {0.0}
}

// Schwarschild time dilation
fn alpha_g(r: &Vector2D) -> f64 {
    sqrt(sqr_alpha_g(r))
}

// Square speed of light (modified by gravitational time dilation)
fn sqr_cprime(r: &Vector2D) -> f64 {
    C2 * sqr_alpha_g(r)
}

// Speed of light (modified by gravitational time dilation)
fn cprime(r: &Vector2D) -> f64 {
    C * alpha_g(r)
}

// Square ellipse metric, the divisor in computing the radius of the
// lightspeed ellipse
fn sqr_ellipse_metric(r: &Vector2D, v: &Vector2D) -> f64 {
    Vector2D::sqr_sin(r, v) * sqr_alpha_g(r) + Vector2D::sqr_cos(r, v)
}

// This represents massive particles along with photons, but in the case of
// photons we set the k2 field (k^2) to f64::MAX as a sentinel (in the case of
// photons we will use refraction formulas to simulate the photon's evolution).
// (As it happens, as a particle's k gets very large the particle will behave
// more light a photon, therefore f64::MAX is a good choice of sentinel.)
//
// These particles will be in free fall relative to the massive body.
pub struct Particle {
    pub r: Vector2D,    // distance to the gravitational body
    pub v: Vector2D,    // velocity relative to the gravitational body
    pub k2: f64         // k^2
                        // "k" factor is the ratio between the particle's total
                        // energy (Pythagorean sum of kinetic and potential)
                        // and the energy needed to escape the gravitational
                        // body, so logically, k = 1 represents escape energy.
                        // The planets have k < 1 relative to the sun; Mercury
                        // has less than Venus, Venus less than Earth, etc.
                        // Voyager 1 has k > 1 (but all of these examples have
                        // k very close to 1).
}

// A simulation should:
//
//  - create a gravitational body (set M above)
//  - create one or more particles or photons with the below constructors
//  - then loop over calling the update() method, recording and/or displaying
//    the particles' positions and velocities as desired (but see also the
//    correct_velocity() method0.
impl Particle {
    // Constructors

    /*
    // XXX DISABLED isotropic physics code
    // Place a particle with zero initial velocity
    pub fn new_vec_place(pos: Vector2D) -> Particle {
        let argk2 = sqr_alpha_g(&pos);

        Particle { r: pos, v: Vector2D::new(0.0, 0.0), k2: argk2 }
    }
    */

    /*
    // XXX DISABLED isotropic physics code
    // Same, but with coordinate arguments instead of Vector2D
    pub fn new_coord_place(xpos: f64, ypos: f64) -> Particle {
        let pos = Vector2D::new(xpos, ypos);
        let argk2 = sqr_alpha_g(&pos);

        Particle { r: pos, v: Vector2D::new(0.0, 0.0), k2: argk2 }
    }
    */

    // Place a particle with given initial velocity
    pub fn new_vec_launch(pos: Vector2D, mut vel: Vector2D) -> Particle {
        // XXX DISABLED isotropic physics code
        /*
        if vel.sqr_magnitude() == 0.0 {
            return Self::new_vec_place(pos);
        }

        let argk2;

        // Compute k2 from position and velocity relative to the massive body
        if sqr_cprime(&pos) - vel.sqr_magnitude() <= 0.0 {
            // Clamp velocity to the local maximum speed of light c', and mark
            // this as a photon
            vel *= sqrt(sqr_cprime(&pos) / vel.sqr_magnitude());
            argk2 = f64::MAX; // it's a photon
        } else {
            // not a photon
            argk2 = sqr_alpha_g(&pos) * sqr_cprime(&pos) / (sqr_cprime(&pos) - vel.sqr_magnitude());
        }

        Particle { r: pos, v: vel, k2: argk2 }
        */

        // If starting velocity is zero, we need to figure out which direction
        // the photon should travel
        if vel.sqr_magnitude() == 0.0 {
            // If at the center, doesn't really matter, leave everything as zero
            if pos.sqr_magnitude() == 0.0 {
                return Particle { r: pos, v: vel, k2: f64::MAX}
            }
            // Otherwise, point the velocity vector at the gravitational body
            else {
                vel.x = -pos.x;
                vel.y = -pos.y;
            }
        }

        // Set the velocity to lightspeed
        if sqr_alpha_g(&pos) == 0.0 {
            vel *= 0.0;
        } else {
            vel *= sqrt(sqr_cprime(&pos) * sqr_alpha_g(&pos) / (sqr_ellipse_metric(&pos, &vel) * vel.sqr_magnitude()));
        }
        Particle { r: pos, v: vel, k2: f64::MAX }
    }

    // Same, but with coordinate arguments instead of Vector2D
    pub fn new_coord_launch(xpos: f64, ypos: f64, vx: f64, vy: f64) -> Particle {
        let pos = Vector2D::new(xpos, ypos);
        let mut vel = Vector2D::new(vx, vy);

        // XXX DISABLED isotropic physics code
        /*
        if vel.sqr_magnitude() == 0.0 {
            return Self::new_vec_place(pos);
        }

        let argk2;

        // See new_vec_launch()
        if sqr_cprime(&pos) - vel.sqr_magnitude() <= 0.0 {
            vel *= sqrt(sqr_cprime(&pos) / vel.sqr_magnitude());
            argk2 = f64::MAX;
        } else {
            argk2 = sqr_alpha_g(&pos) * sqr_cprime(&pos) / (sqr_cprime(&pos) - vel.sqr_magnitude());
        }

        Particle { r: pos, v: vel, k2: argk2 }
        */

        // If starting velocity is zero, we need to figure out which direction
        // the photon should travel
        if vel.sqr_magnitude() == 0.0 {
            // If at the center, doesn't really matter, leave everything as zero
            if pos.sqr_magnitude() == 0.0 {
                return Particle { r: pos, v: vel, k2: f64::MAX}
            }
            // Otherwise, point the velocity vector at the gravitational body
            else {
                vel.x = -pos.x;
                vel.y = -pos.y;
            }
        }

        // Set the velocity to lightspeed
        if sqr_alpha_g(&pos) == 0.0 {
            vel *= 0.0;
        } else {
            vel *= sqrt(sqr_cprime(&pos) * sqr_alpha_g(&pos) / (sqr_ellipse_metric(&pos, &vel) * vel.sqr_magnitude()));
        }
        Particle { r: pos, v: vel, k2: f64::MAX }
    }

    /*
    // XXX DISABLED isotropic physics code
    // In this constructor only the unit vector of velocity is used, and the
    // k^2 argument is used to compute and apply the correct magnitude of the
    // particle's actual velocity vector
    pub fn new_vec_energy(pos: Vector2D, mut vel: Vector2D, argk2: f64) -> Particle {
        if argk2 <= sqr_alpha_g(&pos) || vel.sqr_magnitude() == 0.0 {
            // In the first case, k is less than the floor of k for this
            // location in space, so we ignore k and v and refer to the basic
            // constructor.
            //
            // In the second case, we dont have a direction for velocity, so we
            // ignore the given k and again refer to the basic constructor.
            return Self::new_vec_place(pos);
        }

        // Set the magnitude of velocity based on position and k^2
        vel *= sqrt(sqr_cprime(&pos) * (1.0 - sqr_alpha_g(&pos) / argk2) / vel.sqr_magnitude());

        Particle { r: pos, v: vel, k2: argk2 }
    }
    */

    /*
    // XXX DISABLED isotropic physics code
    // Same, but with coordinate arguments instead of Vector2D
    pub fn new_coord_energy(xpos: f64, ypos: f64, vx: f64, vy: f64, argk2: f64) -> Particle {
        let pos = Vector2D::new(xpos, ypos);
        let mut vel = Vector2D::new(vx, vy);

        if argk2 <= sqr_alpha_g(&pos) || vel.sqr_magnitude() == 0.0 {
            return Self::new_vec_place(pos);
        }

        vel *= sqrt(sqr_cprime(&pos) * (1.0 - sqr_alpha_g(&pos) / argk2) / vel.sqr_magnitude());

        Particle { r: pos, v: vel, k2: argk2 }
    }
    */

    // Single-step (by a small dt) the simulation of the self particle relative
    // to the massive body
    pub fn update(&mut self, dt: f64) {
        if self.r.sqr_magnitude() <= 4.0 * GM * GM / (C2 * C2) {
            // Particles inside the Schwarschild radius do not move (in this
            // simulation).  They are as-if frozen in place.
            self.v.x = 0.0;
            self.v.y = 0.0;
            return;
        }

        // Compute the new position based on current velocity and the time step
        self.r += self.v * dt;

        // A note about terminology:
        //
        //  - "radial" refers to the component of velocity that is parallel to
        //    the vector towards (or away from) the massive body (the radial
        //    vector)
        //
        //  - "tangential" refers to the component of velocity that is
        //    perpendicular to the radial vector
        //
        //  - "axial" refers to the component of acceleration that is parallel
        //    to the velocity vector
        //
        //  - "lateral" refers to the component of acceleration that is
        //    perpendicular to the velocity vector

        // Next we'll compute the instantaneous acceleration due to the massive
        // body and then compute the dv with which to update the particle's
        // velocity.
        //
        // Our theory uses Newtonian gravity, -GM/r^2, modified using special
        // relativistic mechanics.  The axial acceleration for light then turns
        // out to be GM/r^2 (i.e., repulsive!) but a three-term function of r
        // for massive particles.
        //
        // The lateral acceleration is always -GM/r^2, for both light and
        // massive particles (since their _speed_ is zero relative to the
        // massive body if they orbit it in perfect circles, therefore the
        // Newtonian acceleration is not modified by relativistic effects).
        // 
        // The formulas here are -or will be- described in an associated paper.
        //
        // In a variant of this theory the lateral acceleration of light can be
        // of larger magnitude, but never mind that here.

        // This is the magnitude of acceleration along the direction
        // perpendicular to the particle's velocity vector, and it happens to
        // be the Newtonian -GM/r^2.  For massive particles this acceleration
        // is not modified by relativity because their velocity towards or away
        // from the massive body is zero if the particles are moving in
        // circular orbits.  For photons we derived this from refraction with
        // gravitational time dilation (alpha_g) as the slowdown factor.
        // XXX DISABLED isotropic lateral acceleration
        /*
        let lateral_accel = -GM / self.r.sqr_magnitude();
        */
        let lateral_accel = sqr_alpha_g(&self.r) * (Vector2D::sqr_cos(&self.r, &self.v) + sqr_ellipse_metric(&self.r, &self.v)) * GM / (self.r.sqr_magnitude() * pow(sqr_ellipse_metric(&self.r, &self.v), 2.0));

        // XXX DISABLED massive particle case
        /*
        if self.v.sqr_magnitude() == 0.0 {
            // The particle is standing still and will fall directly towards
            // the massive body, and we don't need to consider two acceleration
            // vectors.
            //
            // And note that implicitly self is not a photon here since photons
            // never stand still (except at the Schwarschild radius, which is a
            // case we already handled).
            self.v = self.r * lateral_accel * dt / self.r.magnitude();
            return;
        }   // else
        */

        // XXX DISABLED axial acceleration, use velocity correction at each step
        /*
        let axial_accel;
        if self.k2 != f64::MAX {
            // self is a _massive particle_, *not* a photon, and its
            // acceleration is the Newtonian -GM/r^2 (already computed above)
            // times a relativistic factor which here we'll compute based on
            // self's k, but could be computed from velocity instead, if you
            // want a more complicated formula.
            axial_accel = -lateral_accel * (1.0 - 2.0 / self.k2 + 4.0 * GM / (self.k2 * C2 * self.r.magnitude()));
        } else {
            // else self is a photon.
            //
            // The axial acceleration that light experiences as it moves
            // directly towards to or away from the massive body is GM/r^2
            // (i.e., repulsive) thanks to the slowdown from gravity (alpha_g).
            axial_accel = -lateral_accel;
        }
        */

        // XXX DISABLED cosine and sine computation, now use Vector2D functions
        // Compute the radial and tangential velocity unit vectors at the
        // particle's position (dot and cross product of radius and velocity
        // vectors, respectively, divided by the magnitudes of both vectors).
        // Yields the cosine and sine of the angle between the radius and
        // velocity vectors, respectively. It's presumably faster than using
        // trig functions.
        /*
        let cos = (self.r.x * self.v.x + self.r.y * self.v.y) / sqrt(self.r.sqr_magnitude() * self.v.sqr_magnitude());
        let sin = (self.r.x * self.v.y - self.r.y * self.v.x) / sqrt(self.r.sqr_magnitude() * self.v.sqr_magnitude());
        */

        // XXX DISABLED axial acceleration, use velocity correction at each step
        // Finally, apply the instantaneous accleration to the particle's
        // velocity
        self.v += self.v.norm() * Vector2D::sin(&self.r, &self.v) * lateral_accel * dt / self.v.magnitude();
    }

    // Due to the time steps not being infinitessimal we accumulate some
    // error in the particle's velocity in each call to update().  This method
    // corrects for that.  We do not invoke it in the update() method because
    // we are not certain whether it is computationally more efficient to
    // correct the velocity at each time step or to simply use smaller time
    // steps. Due to the possibility of floating point errors when the time
    // steps get too small, this correction function may only be efficient when
    // the time steps are already very small.
    pub fn correct_velocity(&mut self) {
        if self.k2 <= sqr_alpha_g(&self.r) || self.v.sqr_magnitude() == 0.0 {
            return;
        }

        // XXX DISABLED unified velocity correction for both massive and massless particles
        // Set the magnitude of velocity based on position and k^2
        /*
        self.v *= sqrt(sqr_cprime(&self.r) * (1.0 - sqr_alpha_g(&self.r) / self.k2) / self.v.sqr_magnitude());
        */
        self.v *= sqrt(sqr_cprime(&self.r) * sqr_alpha_g(&self.r) / (sqr_ellipse_metric(&self.r, &self.v) * self.v.sqr_magnitude()));
    }
}
