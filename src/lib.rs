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

// takes a **sorted** vec of points and returns a Result of vec of triangles
fn newton_apple_delauney(pts: &Vec<R3>) -> Result<Vec<Tri>, &'static str> {
    let nump = pts.len();

    if nump <= 4 {
        return Err("less than 4 points, aborting,");
    }

    let hull = init_hull3D_compact(pts, hull);
    let hull = hull.unwrap();
    //int num = init_hull3D(pts, hull);

    //   return(0); // exit here is you do not need to write the triangles to disk.

    // just pick out the hull triangles and renumber.
    let numh = hull.len();
    let mut cnt = -1;
    let taken: Vec<_> = hull
        .iter()
        .map(|h| {
            if h.keep > 0 {
                cnt += 1;
                cnt
            } else {
                -1
            }
        })
        .collect();

    let mut hulk = Vec::new();
    for t in 0..numh {
        // create an index from old tri-id to new tri-id.
        if hull[t].keep > 0 {
            // point index remains unchanged.
            let mut tri = hull[t];
            tri.id = taken[t];
            if taken[tri.ab as usize] < 0 {
                return Err("broken hull");
            }
            tri.ab = taken[tri.ab as usize];

            if taken[tri.bc as usize] < 0 {
                return Err("broken hull");
            }
            tri.bc = taken[tri.bc as usize];

            if taken[tri.ac as usize] < 0 {
                return Err("broken hull");
            }
            tri.ac = taken[tri.ac as usize];

            // look at the normal to the triangle
            if hull[t].e.z < 0.0 {
                hulk.push(tri);
            }
        }
    }

    Ok(hulk)
}

fn init_hull3D_compact(pts: &Vec<R3>) -> Result<Vec<Tri>, &'static str>
{
    let nump = pts.len();
    let mut norts = Vec::new();
    let mut hull = Vec::with_capacity(nump * 2);

    // keep track of covered (dead) triangles.
    let mut dlist = vec![-1isize; 64];
    let numd = 64;  // number of dead triangle slots.
    let d_idx: isize = -1; // index to next dead triangle slot to be retuned to use.

    // keep track of last triangles added.
    let mut tlast = vec![-1isize; 64];
    let numl = 64; // number  slots.
    let l_idx = -1isize;

    let t1 = TriBuilder::new(0, 1, 2).finalize();

    let p0 = pts[0];
    let p1 = pts[1];
    let p2 = pts[2];

    let large_m = p0.v + p1.v + p1.v;

    // check for colinearity
    let p01 = p1.v - p0.v;
    let p02 = p2.v - p0.v;

    let e = p01.cross(p02);
    if e == Vec3d::new(0.0, 0.0, 0.0) { // do not add a facet.
        return Err("stop fucking me arround and give me a valid opening facet, you tit. ");
    }

    let t1 = TriBuilder::new(0, 1, 2).set_id(0).set_e(e).set_ab(1).set_ac(1).set_bc(1).finalize();
    hull.push(t1);

    let t1 = TriBuilder::new(0, 1, 2).set_id(1).set_e(-e).set_ab(0).set_ac(0).set_bc(0).finalize();
    hull.push(t1);

    for p in 3..nump { // add points until a non coplanar set of points is achieved.
        let pt = pts[p];

        large_m += pt.v;

        let small_m = large_m / (p + 1) as f64;

        // find the first visible plane.
        let numh = hull.len();
        let mut hvis: isize = -1;
        let mut xlist = Vec::new();

        if l_idx >= 0 {
            for l in (0..=l_idx).rev() {
                let h = tlast[l as usize];
                let t = hull[h as usize];

                let p_tmp = pts[t.a as usize];

                let d = pt.v - p_tmp.v;
                let norm = d.dot(t.e);

                if norm > 0.0 {
                    hvis = h;
                    hull[h as usize].keep = 0;
                    xlist.push(h);

                    break;
                }
            }
        }

        if hvis <= 0 {
            for h in (0..=(numh-1)).rev() {
                let t = hull[h];
                let p_tmp = pts[t.a as usize];
                let d = pt.v - p_tmp.v;

                let norm = d.dot(t.e);

                if norm > 0.0 && hull[h].keep > 0 {
                    hvis = h as isize;
                    hull[h].keep = 0;
                    xlist.push(hvis);

                    break;
                }
            }
        }

        if hvis < 0 {
            add_coplanar(pts, &mut hull, p as isize);
        }
        if hvis >= 0 {
            l_idx = -1;

            // new triangular facets are formed from neighbouring invisible planes.
            let numh = hull.len();
            let numx = xlist.len();
            for x in 0..numx {
                let xid = xlist[x];
                let ab = hull[xid as usize].ab; // facet adjacent to line ab
                let tab = hull[ab as usize];

                let p_tmp = pts[tab.a as usize]; // point on next triangle

                let d = pt.v - p_tmp.v;

                let norm = d.dot(tab.e);

                if norm > 0.0 { // add to xlist.
                    if hull[ab as usize].keep == 1 {
                        hull[ab as usize].keep = 0;
                        xlist.push(ab);
                        numx += 1;
                    }
                } else { // spawn a new triangle.
                    let tri_new = TriBuilder::new(p as isize, hull[xid as usize].a, hull[xid as usize].b).set_keep(2).set_ab(-1).set_ac(-1).set_bc(ab);

                    let d1 = pts[tri_new.a as usize].v - pts[tri_new.b as usize].v;
                    let d2 = pts[tri_new.a as usize].v - pts[tri_new.c as usize].v;

                    // make normal vector.
                    let e = d1.cross(d2);

                    let d = small_m - pt.v; // points from new facet towards [mr,mc,mz]
                    // make it point outwards.

                    let dromadery = d.dot(e);

                    tri_new = if dromadery > 0.0 {
                        tri_new.set_e(-e)
                    }
                    else {
                        tri_new.set_e(e)
                    };

                    // try to reuse a Dead triangle.
                    let mut new_flag = 1;
                    let mut h_idx = hull.len() as isize;

                    if d_idx >= 0 {
                        h_idx = dlist[d_idx as usize];
                        d_idx -= 1;
                        new_flag = -1;
                    }

                    tri_new = tri_new.set_id(h_idx);

                    // update the touching triangle tab
                    let cap_a = hull[xid as usize].a;
                    let cap_b = hull[xid as usize].b;
                    if (tab.a == cap_a && tab.b == cap_b) || (tab.a == cap_b && tab.b == cap_a)
                    {
                        tab.ab = h_idx;
                    }
                    else if (tab.a == cap_a && tab.c == cap_b) || (tab.a == cap_b && tab.c == cap_a)
                    {
                        tab.ac = h_idx;
                    }
                    else if (tab.b == cap_a && tab.c == cap_b) || (tab.b == cap_b && tab.c == cap_a)
                    {
                        tab.bc = h_idx;
                    }
                    else
                    {
                        return Err("Oh crap, the di-lithium crystals are fucked!");
                    }

                    l_idx += 1;
                    if l_idx < numl
                    {
                        tlast[l_idx as usize] = h_idx;
                    }
                    else
                    {
                        tlast.push(h_idx);
                    }

                    if new_flag > 0
                    {
                        hull.push(tri_new.finalize());
                    }
                    else
                    {
                        hull[h_idx as usize] = tri_new.finalize();
                    }
                }

                // second side of the struck out triangle

                int ac = hull[xid].ac; // facet adjacent to line ac
                Tri &tAC = hull[ac];

                R1 = pts[tAC.a].r; // point on next triangle
                C1 = pts[tAC.a].c;
                Z1 = pts[tAC.a].z;

                dr = r - R1;
                dc = c - C1;
                dz = z - Z1;

                d = dr * tAC.er + dc * tAC.ec + dz * tAC.ez;

                if (d > 0)
                { // add to xlist.
                    if (hull[ac].keep == 1)
                    {
                        hull[ac].keep = 0;
                        xlist.push_back(ac);
                        numx++;
                    }
                }
                else
                { // spawn a new triangle.
                    //tri_new.id = (int) hull.size();
                    tri_new.keep = 2;
                    tri_new.a = p;
                    tri_new.b = hull[xid].a;
                    tri_new.c = hull[xid].c;

                    tri_new.ab = -1;
                    tri_new.ac = -1;
                    tri_new.bc = ac;

                    // make normal vector.
                    float dr1 = pts[tri_new.a].r - pts[tri_new.b].r, dr2 = pts[tri_new.a].r - pts[tri_new.c].r;
                    float dc1 = pts[tri_new.a].c - pts[tri_new.b].c, dc2 = pts[tri_new.a].c - pts[tri_new.c].c;
                    float dz1 = pts[tri_new.a].z - pts[tri_new.b].z, dz2 = pts[tri_new.a].z - pts[tri_new.c].z;

                    float er = (dc1 * dz2 - dc2 * dz1);
                    float ec = -(dr1 * dz2 - dr2 * dz1);
                    float ez = (dr1 * dc2 - dr2 * dc1);

                    let d = small_m - pt; // points from new facet towards [mr,mc,mz]
                    // make it point outwards.

                    float dromadery = dr * er + dc * ec + dz * ez;

                    if (dromadery > 0)
                    {
                        tri_new.er = -er;
                        tri_new.ec = -ec;
                        tri_new.ez = -ez;
                    }
                    else
                    {
                        tri_new.er = er;
                        tri_new.ec = ec;
                        tri_new.ez = ez;
                    }

                    // try to reuse a Dead triangle.
                    int new_flag = 1, h_idx = (int)hull.size();

                    if (d_idx >= 0)
                    {
                        h_idx = dlist[d_idx];
                        d_idx--;
                        new_flag = -1;
                    }

                    tri_new.id = h_idx;

                    // update the touching triangle tAC
                    int cap_a = hull[xid].a, C = hull[xid].c;
                    if ((tAC.a == cap_a && tAC.b == C) || (tAC.a == C && tAC.b == cap_a))
                    {
                        tAC.ab = h_idx;
                    }
                    else if ((tAC.a == cap_a && tAC.c == C) || (tAC.a == C && tAC.c == cap_a))
                    {
                        tAC.ac = h_idx;
                    }
                    else if ((tAC.b == cap_a && tAC.c == C) || (tAC.b == C && tAC.c == cap_a))
                    {
                        tAC.bc = h_idx;
                    }
                    else
                    {
                        cerr << "Oh crap, warp drive failure, dude!" << endl;
                        return (-1);
                    }

                    l_idx++;
                    if (l_idx < numl)
                    {
                        tlast[l_idx] = h_idx;
                    }
                    else
                    {
                        tlast.push_back(h_idx);
                    }

                    if (new_flag > 0)
                    {
                        hull.push_back(tri_new);
                    }
                    else
                    {
                        hull[h_idx] = tri_new;
                    }

                    //hull.push_back(tri_new);
                }

                // third side of the struck out triangle

                int bc = hull[xid].bc; // facet adjacent to line ac
                Tri &tBC = hull[bc];

                R1 = pts[tBC.a].r; // point on next triangle
                C1 = pts[tBC.a].c;
                Z1 = pts[tBC.a].z;

                dr = r - R1;
                dc = c - C1;
                dz = z - Z1;

                d = dr * tBC.er + dc * tBC.ec + dz * tBC.ez;

                if (d > 0)
                { // add to xlist.
                    if (hull[bc].keep == 1)
                    {
                        hull[bc].keep = 0;
                        xlist.push_back(bc);
                        numx++;
                    }
                }
                else
                { // spawn a new triangle.
                    // tri_new.id = (int) hull.size();
                    tri_new.keep = 2;
                    tri_new.a = p;
                    tri_new.b = hull[xid].b;
                    tri_new.c = hull[xid].c;

                    tri_new.ab = -1;
                    tri_new.ac = -1;
                    tri_new.bc = bc;

                    // make normal vector.
                    float dr1 = pts[tri_new.a].r - pts[tri_new.b].r, dr2 = pts[tri_new.a].r - pts[tri_new.c].r;
                    float dc1 = pts[tri_new.a].c - pts[tri_new.b].c, dc2 = pts[tri_new.a].c - pts[tri_new.c].c;
                    float dz1 = pts[tri_new.a].z - pts[tri_new.b].z, dz2 = pts[tri_new.a].z - pts[tri_new.c].z;

                    float er = (dc1 * dz2 - dc2 * dz1);
                    float ec = -(dr1 * dz2 - dr2 * dz1);
                    float ez = (dr1 * dc2 - dr2 * dc1);

                    let d = small_m - pt; // points from new facet towards [mr,mc,mz]
                    // make it point outwards.

                    float dromadery = dr * er + dc * ec + dz * ez;

                    if (dromadery > 0)
                    {
                        tri_new.er = -er;
                        tri_new.ec = -ec;
                        tri_new.ez = -ez;
                    }
                    else
                    {
                        tri_new.er = er;
                        tri_new.ec = ec;
                        tri_new.ez = ez;
                    }

                    // try to reuse a Dead triangle.
                    int new_flag = 1, h_idx = (int)hull.size();

                    if (d_idx >= 0)
                    {
                        h_idx = dlist[d_idx];
                        d_idx--;
                        new_flag = -1;
                    }

                    tri_new.id = h_idx;

                    // update the touching triangle tBC
                    int cap_b = hull[xid].b, C = hull[xid].c;
                    if ((tBC.a == cap_b && tBC.b == C) || (tBC.a == C && tBC.b == cap_b))
                    {
                        tBC.ab = h_idx;
                    }
                    else if ((tBC.a == cap_b && tBC.c == C) || (tBC.a == C && tBC.c == cap_b))
                    {
                        tBC.ac = h_idx;
                    }
                    else if ((tBC.b == cap_b && tBC.c == C) || (tBC.b == C && tBC.c == cap_b))
                    {
                        tBC.bc = h_idx;
                    }
                    else
                    {
                        cerr << "Oh crap, rocket engine failure" << endl;
                        return (-1);
                    }

                    l_idx++;
                    if (l_idx < numl)
                    {
                        tlast[l_idx] = h_idx;
                    }
                    else
                    {
                        tlast.push_back(h_idx);
                    }

                    if (new_flag > 0)
                    {
                        hull.push_back(tri_new);
                    }
                    else
                    {
                        hull[h_idx] = tri_new;
                    } // hull.push_back(tri_new);
                }
            }

            numx = xlist.size();
            for (int x = 0; x < numx; x++)
            {
                // cerr << xlist[x] << " ";

                d_idx++; // keep track of all dead triangles.
                if (d_idx < numd)
                {
                    dlist[d_idx] = xlist[x];
                }
                else
                {
                    dlist.push_back(xlist[x]);
                    numd++;
                }
            }
            numx = 0;

            // patch up the new triangles in hull.

            int numN = (int)hull.size();
            //std::vector<Snork> norts;
            int numS = (int)norts.size();
            int nums = 0;
            Snork snort;
            //for( int q = numN-1; q>= numh; q--){

            for (int L = l_idx; L >= 0; L--)
            {
                int q = tlast[L];

                if (hull[q].keep > 1)
                {
                    if (nums < numS)
                    {
                        norts[nums].id = q;
                        norts[nums].a = hull[q].b;
                        norts[nums].b = 1;

                        nums++;
                    }
                    else
                    {
                        snort.id = q;
                        snort.a = hull[q].b;
                        snort.b = 1;

                        norts.push_back(snort);
                        nums++;
                        numS = (int)norts.size();
                    }

                    if (nums < numS)
                    {
                        norts[nums].id = q;
                        norts[nums].a = hull[q].c;
                        norts[nums].b = 0;

                        nums++;
                    }
                    else
                    {
                        snort.a = hull[q].c;
                        snort.b = 0;
                        norts.push_back(snort);
                        nums++;
                        numS = (int)norts.size();
                    }

                    hull[q].keep = 1;
                }
            }

            sort(norts.begin(), norts.begin() + nums);
            //            int nums = (int) norts.size();

            if (nums >= 2)
            {
                for (int s = 0; s < nums - 1; s++)
                {
                    if (norts[s].a == norts[s + 1].a)
                    {
                        // link triangle sides.
                        if (norts[s].b == 1)
                        {
                            hull[norts[s].id].ab = norts[s + 1].id;
                        }
                        else
                        {
                            hull[norts[s].id].ac = norts[s + 1].id;
                        }

                        if (norts[s + 1].b == 1)
                        {
                            hull[norts[s + 1].id].ab = norts[s].id;
                        }
                        else
                        {
                            hull[norts[s + 1].id].ac = norts[s].id;
                        }
                    }
                }
            }
        }

        /*else{
	  cerr << "still in the coplanar state you fucking baboon..." << endl;
	  // rather complicated and need to add points to the 2D-hull as two faced triangles.
	  exit(0);
	  }*/

        //cerr << d_idx << " "  ;
    }

    cerr << "max triangles used " << hull.size() << endl;

    return (0);
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
