use std::cmp::Ordering;
use std::ops::{Add, Mul, Sub, Neg};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3d {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3d {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3d { x, y, z }
    }

    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    pub fn norm_sqr(self) -> f64 {
        self.dot(self)
    }

    pub fn norm(self) -> f64 {
        f64::sqrt(self.norm_sqr())
    }

    pub fn distance(self, rhs: Self) -> f64 {
        (self - rhs).norm()
    }

    pub fn cross(self, rhs: Self) -> Self {
        Vec3d::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x,
        )
    }
}

impl Add for Vec3d {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vec3d {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f64> for Vec3d {
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Neg for Vec3d {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Mul<Vec3d> for f64 {
    type Output = Vec3d;

    fn mul(self, other: Vec3d) -> Self::Output {
        other.mul(self)
    }
}

// // sort into descending order (for use in corner responce ranking).
// inline bool operator<(const R3 &a, const R3 &b)
// {
//   if (a.z == b.z)
//   {
//     if (a.r == b.r)
//     {
//       return a.c < b.c;
//     }
//     return a.r < b.r;
//   }
//   return a.z < b.z;
// };

impl PartialOrd for Vec3d {
    fn partial_cmp(&self, other: &Vec3d) -> Option<Ordering> {
        if self.z == other.z {
            if self.y == other.y {
                if self.x == other.x {
                    return Some(Ordering::Equal);
                } else if self.x < other.x {
                    return Some(Ordering::Less);
                } else if self.x > other.x {
                    return Some(Ordering::Greater);
                } else {
                    return None;
                }
            }
            if self.y < other.y {
                return Some(Ordering::Less);
            } else if self.y > other.y {
                return Some(Ordering::Greater);
            } else {
                return None;
            }
        }
        if self.z < other.z {
            return Some(Ordering::Less);
        } else if self.z > other.z {
            return Some(Ordering::Greater);
        } else {
            return None;
        }
    }
}

pub fn sort_and_deduplicate(mut points: Vec<Vec3d>) -> Vec<Vec3d> {
    points.sort_by(|a, b| a.partial_cmp(b).unwrap());
    points.dedup();
    points
}

pub fn cross_test(
    points: &Vec<Vec3d>,
    p1: usize,
    p2: usize,
    p3: usize,
    id: usize,
) -> (isize, Vec3d) {
    let a = points[p1];
    let b = points[p2];
    let c = points[p3];
    let x = points[id];

    let ab = b - a;
    let ac = b - a;
    let ax = b - a;

    let e = ab.cross(ax);
    let k = ab.cross(ac);

    let globit = k.dot(e);

    if globit > 0.0 {
        return (1, e);
    } else if globit == 0.0 {
        return (0, e);
    } else {
        return (-1, e);
    }
}
