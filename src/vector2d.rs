use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Div;
use std::ops::DivAssign;
use libm::sqrt;
use libm::pow;

#[derive(Copy, Clone)]
pub struct Vector2D {pub x: f64, pub y: f64}

impl Vector2D {
    pub fn zero() -> Self {
        Self {x: 0.0, y: 0.0}
    }

    pub fn new(xval: f64, yval: f64) -> Self {
        Self {x: xval, y: yval}
    }

    pub fn norm(&self) -> Self {
        Self {x: self.y, y: -self.x}
    }

    pub fn sqr_magnitude(&self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    pub fn magnitude(&self) -> f64 {
        sqrt(self.sqr_magnitude())
    }

    pub fn sqr_sin(a: &Self, b: &Self) -> f64 {
        pow(a.x * b.y - a.y * b.x, 2.0) / (a.sqr_magnitude() * b.sqr_magnitude())
    }

    pub fn sin(a: &Self, b: &Self) -> f64 {
        sqrt(Self::sqr_sin(a, b))
    }

    pub fn sqr_cos(a: &Self, b: &Self) -> f64 {
        pow(a.x * b.x + a.y * b.y, 2.0) / (a.sqr_magnitude() * b.sqr_magnitude())
    }

    pub fn cos(a: &Self, b: &Self) -> f64 {
        sqrt(Self::sqr_cos(a, b))
    }
}

impl Mul<f64> for Vector2D {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs)
    }
}

impl MulAssign<f64> for Vector2D {
    fn mul_assign(&mut self, rhs: f64) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

impl Div<f64> for Vector2D {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self::new(self.x / rhs, self.y / rhs)
    }
}

impl DivAssign<f64> for Vector2D {
    fn div_assign(&mut self, rhs: f64) {
        self.x /= rhs;
        self.y /= rhs;
    }
}

impl Add for Vector2D {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl AddAssign for Vector2D {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}
