mod snork;
mod vec3d;

use snork::Snork;
use std::mem;
pub use vec3d::{sort_and_deduplicate, Vec3d};

struct TriBuilder {
    id: isize,
    keep: isize,
    a: isize,
    b: isize,
    c: isize,
    ab: isize,
    bc: isize,
    ac: isize,
    e: Vec3d,
}

impl TriBuilder {
    fn new(a: isize, b: isize, c: isize) -> Self {
        TriBuilder {
            id: 0,
            keep: 1,
            a,
            b,
            c,
            ab: -1,
            bc: -1,
            ac: -1,
            e: Vec3d::new(0.0, 0.0, 0.0),
        }
    }

    fn set_id(mut self, id: isize) -> Self {
        self.id = id;
        self
    }

    fn set_keep(mut self, keep: isize) -> Self {
        self.keep = keep;
        self
    }

    fn set_ab(mut self, ab: isize) -> Self {
        self.ab = ab;
        self
    }

    fn set_ac(mut self, ac: isize) -> Self {
        self.ac = ac;
        self
    }

    fn set_bc(mut self, bc: isize) -> Self {
        self.bc = bc;
        self
    }

    fn set_e(mut self, e: Vec3d) -> Self {
        self.e = e;
        self
    }

    fn finalize(self) -> Tri {
        Tri {
            id: self.id,
            keep: self.keep,
            a: self.a,
            b: self.b,
            c: self.c,
            ab: self.ab,
            bc: self.bc,
            ac: self.ac,
            e: self.e,
        }
    }
}

struct Tri {
    id: isize,
    keep: isize,
    a: isize,
    b: isize,
    c: isize,
    ab: isize,
    bc: isize,
    ac: isize,
    e: Vec3d,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct R3 {
    v: Vec3d,
    id: usize,
}

impl R3 {
    pub fn new(v: Vec3d, id: usize) -> Self {
        R3 { v, id }
    }

    pub fn dot(self, rhs: Self) -> f64 {
        self.v.dot(rhs.v)
    }

    pub fn norm_sqr(self) -> f64 {
        self.dot(self)
    }

    pub fn norm(self) -> f64 {
        f64::sqrt(self.norm_sqr())
    }

    pub fn distance(self, rhs: Self) -> f64 {
        (self.v - rhs.v).norm()
    }

    pub fn cross(self, rhs: Self) -> Vec3d {
        self.v.cross(rhs.v)
    }
}

pub fn from_vec3d_to_r3(points: Vec<Vec3d>) -> Vec<R3> {
    points
        .into_iter()
        .enumerate()
        .map(|(i, p)| R3::new(p, i))
        .collect()
}

pub fn from_ref_r3_to_vec3d(points: &Vec<R3>) -> Vec<Vec3d> {
    points.iter().map(|r| r.v).collect::<Vec<Vec3d>>()
}

fn add_coplanar(pts: &Vec<R3>, hull: &mut Vec<Tri>, id: isize) {
    let numh = hull.len();
    for k in 0..numh {
        //find vizible edges. from external edges.
        if hull[k].c == hull[hull[k].ab as usize].c {
            // ->  ab is an external edge.
            // test this edge for visibility from new point pts[id].
            let p1 = hull[k].a;
            let p2 = hull[k].b;
            let p3 = hull[k].c;

            let (zot, e) = vec3d::cross_test(
                &from_ref_r3_to_vec3d(pts),
                p1 as usize,
                p2 as usize,
                p3 as usize,
                id as usize,
            );

            if zot < 0 {
                // visible edge facet, create 2 new hull plates.
                let xx = hull[k].e.dot(e);
                let mut up = TriBuilder::new(id, p1, p2)
                    .set_e(e)
                    .set_id(hull.len() as isize)
                    .set_ab(-1)
                    .set_ac(-1)
                    .set_keep(2);
                let mut down = TriBuilder::new(id, p1, p2)
                    .set_e(-e)
                    .set_id(hull.len() as isize + 1)
                    .set_ab(-1)
                    .set_ac(-1)
                    .set_keep(2);

                if xx > 0.0 {
                    up = up.set_bc(k as isize);
                    down = down.set_bc(hull[k].ab);

                    hull[k].ab = up.id;
                    hull[down.bc as usize].ab = down.id;
                } else {
                    down = down.set_bc(k as isize);
                    up = up.set_bc(hull[k].ab);

                    hull[k].ab = down.id;
                    hull[up.bc as usize].ab = up.id;
                }

                hull.push(up.finalize());
                hull.push(down.finalize());
            }
        }

        if hull[k].a == hull[hull[k].bc as usize].a {
            // bc is an external edge.
            // test this edge for visibility from new point pts[id].
            let p1 = hull[k].b;
            let p2 = hull[k].c;
            let p3 = hull[k].a;

            let (zot, e) = vec3d::cross_test(
                &from_ref_r3_to_vec3d(pts),
                p1 as usize,
                p2 as usize,
                p3 as usize,
                id as usize,
            );

            if zot < 0 {
                // visible edge facet, create 2 new hull plates.
                let mut up = TriBuilder::new(id, p1, p2)
                    .set_id(hull.len() as isize)
                    .set_e(e)
                    .set_keep(2)
                    .set_ab(-1)
                    .set_ac(-1);
                let mut down = TriBuilder::new(id, p1, p2)
                    .set_id(hull.len() as isize + 1)
                    .set_e(-e)
                    .set_keep(2)
                    .set_ab(-1)
                    .set_ac(-1);

                let xx = hull[k].e.dot(e);
                if xx > 0.0 {
                    up = up.set_bc(k as isize);
                    down = down.set_bc(hull[k].bc);

                    hull[k].bc = up.id;
                    hull[down.bc as usize].bc = down.id;
                } else {
                    down = down.set_bc(k as isize);
                    up = up.set_bc(hull[k].bc);

                    hull[k].bc = down.id;
                    hull[up.bc as usize].bc = up.id;
                }

                hull.push(up.finalize());
                hull.push(down.finalize());
            }
        }

        if hull[k].b == hull[hull[k].ac as usize].b {
            // ac is an external edge.
            // test this edge for visibility from new point pts[id].
            let p1 = hull[k].a;
            let p2 = hull[k].c;
            let p3 = hull[k].b;

            let (zot, e) = vec3d::cross_test(
                &from_ref_r3_to_vec3d(pts),
                p1 as usize,
                p2 as usize,
                p3 as usize,
                id as usize,
            );

            if zot < 0 {
                // visible edge facet, create 2 new hull plates.
                let mut up = TriBuilder::new(id, p1, p2)
                    .set_id(hull.len() as isize)
                    .set_e(e)
                    .set_keep(2)
                    .set_ab(-1)
                    .set_ac(-1);
                let mut down = TriBuilder::new(id, p1, p2)
                    .set_id(hull.len() as isize + 1)
                    .set_e(-e)
                    .set_keep(2)
                    .set_ab(-1)
                    .set_ac(-1);

                let xx = hull[k].e.dot(e);
                if xx > 0.0 {
                    up = up.set_bc(k as isize);
                    down = down.set_bc(hull[k].ac);

                    hull[k].ac = up.id;
                    hull[down.bc as usize].ac = down.id;
                } else {
                    down = down.set_bc(k as isize);
                    up = up.set_bc(hull[k].ac);

                    hull[k].ac = down.id;
                    hull[up.bc as usize].ac = up.id;
                }

                hull.push(up.finalize());
                hull.push(down.finalize());
            }
        }
    }

    // fix up the non asigned hull adjecencies (correctly).

    let num_n = hull.len();
    let mut norts = Vec::new();
    for q in (numh..=(num_n - 1)).rev() {
        if hull[q].keep > 1 {
            norts.push(Snork::new(q as isize, hull[q].b, 1));
            norts.push(Snork::new(q as isize, hull[q].c, 0));

            hull[q].keep = 1;
        }
    }

    norts.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let nums = norts.len();
    norts.push(Snork::new(-1, -1, -1));
    norts.push(Snork::new(-2, -2, -2));

    if nums >= 2 {
        let mut s = 0;
        while s < (nums - 1) {
            if norts[s].a == norts[s + 1].a {
                // link triangle sides.
                if norts[s].a != norts[s + 2].a {
                    // edge of figure case
                    if norts[s].b == 1 {
                        hull[norts[s].id as usize].ab = norts[s + 1].id;
                    } else {
                        hull[norts[s].id as usize].ac = norts[s + 1].id;
                    }

                    if norts[s + 1].b == 1 {
                        hull[norts[s + 1].id as usize].ab = norts[s].id;
                    } else {
                        hull[norts[s + 1].id as usize].ac = norts[s].id;
                    }
                    s += 1;
                } else {
                    // internal figure boundary 4 junction case.
                    let mut s1 = s + 1;
                    let mut s2 = s + 2;
                    let mut s3 = s + 3;

                    let id = norts[s].id;
                    let mut id1 = norts[s1].id;
                    let mut id2 = norts[s2].id;
                    let mut id3 = norts[s3].id;

                    // check normal directions of id and id1..3
                    let mut barf = hull[id as usize].e.dot(hull[id1 as usize].e);
                    if barf <= 0.0 {
                        barf = hull[id as usize].e.dot(hull[id2 as usize].e);
                        if barf > 0.0 {
                            mem::swap(&mut id1, &mut id2);
                            mem::swap(&mut s1, &mut s2);
                        } else {
                            barf = hull[id as usize].e.dot(hull[id3 as usize].e);
                            if barf > 0.0 {
                                mem::swap(&mut id1, &mut id3);
                                mem::swap(&mut s1, &mut s3);
                            }
                        }
                    }

                    if norts[s].b == 1 {
                        hull[norts[s].id as usize].ab = norts[s1].id;
                    } else {
                        hull[norts[s].id as usize].ac = norts[s1].id;
                    }

                    if norts[s1].b == 1 {
                        hull[norts[s1].id as usize].ab = norts[s].id;
                    } else {
                        hull[norts[s1].id as usize].ac = norts[s].id;
                    }

                    // use s2 and s3

                    if norts[s2].b == 1 {
                        hull[norts[s2].id as usize].ab = norts[s3].id;
                    } else {
                        hull[norts[s2].id as usize].ac = norts[s3].id;
                    }

                    if norts[s3].b == 1 {
                        hull[norts[s3].id as usize].ab = norts[s2].id;
                    } else {
                        hull[norts[s3].id as usize].ac = norts[s2].id;
                    }

                    s += 3;
                }
            }
        }
    }
}

// // from NewtonApple_hull3D.cpp

// int read_R3(std::vector<R3> &pts, char *fname);
// void write_R3(std::vector<R3> &pts, char *fname);
// void write_Tris(std::vector<Tri> &ts, char *fname);
// void write_Tris_stl(std::vector<Tri> &ts, std::vector<R3> &pts, char *fname);

// int NewtonApple_Delaunay(std::vector<R3> &pts, std::vector<Tri> &hulk);
// int NewtonApple_hull_3D(std::vector<R3> &pts, std::vector<Tri> &hull);

// int init_hull3D(std::vector<R3> &pts, std::vector<Tri> &hull);
// int init_hull3D_compact(std::vector<R3> &pts, std::vector<Tri> &hull);

// #endif
