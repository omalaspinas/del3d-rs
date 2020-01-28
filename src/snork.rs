use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Snork {
    pub id: isize,
    pub a: isize,
    pub b: isize,
}

impl Snork {
    pub fn new(id: isize, a: isize, b: isize) -> Self {
        Snork { id, a, b }
    }
}

impl PartialOrd for Snork {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.a == other.a {
            if self.b == other.b {
                return Some(Ordering::Equal);
            } else if self.b < other.b {
                return Some(Ordering::Less);
            } else if self.b > other.b {
                return Some(Ordering::Greater);
            } else {
                return None;
            }
        }
        if self.a < other.a {
            return Some(Ordering::Less);
        } else if self.a > other.a {
            return Some(Ordering::Greater);
        } else {
            return None;
        }
    }
}
