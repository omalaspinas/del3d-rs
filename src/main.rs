use del3d_rs::{from_vec3d_to_r3};
use del3d_rs::{sort_and_deduplicate, Vec3d};
use rand::prelude::*;

// fn deduplicate( pts: Vec<R3> ) -> (Vec<R3>, Vec<usize>) {

//     let nump = (int) pts.size();
//     std::vector<R3> dpx;
//     R3 d;
//     for( int k=0; k<nump; k++){
//         d.r = pts[k].r;
//         d.c = pts[k].c;
//         d.z = pts[k].z;
//         d.id = k;
//         dpx.push_back(d);
//     }

//     sort(dpx.begin(), dpx.end());

//     cerr << "de-duplicating ";
//     pts2.clear();
//     pts2.push_back(pts[dpx[0].id]);
//     pts2[0].id = 0;
//     int cnt = 1;

//     for( int k=0; k<nump-1; k++){
//         if( dpx[k].r == dpx[k+1].r && dpx[k].c == dpx[k+1].c && dpx[k].z == dpx[k+1].z ){
//             outx.push_back( dpx[k+1].id);
//         }
//         else{
//             pts[dpx[k+1].id].id = cnt;
//             pts2.push_back(pts[dpx[k+1].id]);
//             cnt++;
//         }
//     }

//     cerr << "removed  " << outx.size() << " points " << endl;

//     return((int) outx.size());
// }

// #[derive(Debug, Clone, Copy, PartialEq)]
// struct Tri {
//     id: usize,
//     keep: usize,
//     a: usize,
//     b: usize,
//     c: usize,
//     ab: usize, // adjacent edges index to neighbouring triangle.
//     bc: usize,
//     ac: usize,
//     ex: f64, // visible normal to triangular facet.
//     ey: f64,
//     ez: f64,
// }

// fn main() {

//     let num_points = 100;
//     let mut rng = rand::thread_rng();
//     let points = (0..num_points).into_iter().enumerate().map(|(id, _)| {
//         let x = rng.gen();
//         let y = rng.gen();
//         let z = x*x + y*y;

//         R3::new(x, y, z, id)
//     });

//         std::vector<Tri> tris;

//         std::vector<int> outx;
//         int nx = de_duplicateR3( pts, outx, pts2);
//         pts.clear();

//         write_R3(pts2, "pts.mat");
//         cerr << pts2.size() << " randomly generated points int R2/R3 written to pts.mat" << endl;

//         struct timeval tv1, tv2; // slight swizzle as pt set is now sorted.
//         gettimeofday(&tv1, NULL);

//         int ts = NewtonApple_Delaunay( pts2, tris);
//         // int ts = NewtonApple_hull_3D( pts2, tris);

//         gettimeofday(&tv2, NULL);
//         float tx =  (tv2.tv_sec + tv2.tv_usec / 1000000.0) - ( tv1.tv_sec + tv1.tv_usec / 1000000.0);
//         write_Tris_stl(tris, pts2, "triangles.stl");
//         pts2.clear();

//         cerr <<  tx << " seconds for triangulation" << endl;

//         write_Tris(tris, "triangles.mat");
//         cerr << tris.size() << " triangles written to triangles.mat" << endl;

//         exit(0);
// }

fn main() {
    let num_points = 10;
    let ary = [ -17.9199, 5.34844, -24.185,
                      -17.166, -4.95278, -18.5105,
                      14.922, 20.5824, -15.1224,
                      -19.5596, 24.9462, -14.0872,
                      -8.23886, 13.4115, -11.1113,
                      2.6985, -1.13015, 6.44354,
                      20.8098, 6.78559, 10.8648,
                      17.0094, -5.28085, 14.155,
                      -12.8557, -18.1384, 15.2088,
                      -6.76078, 0.670046, 22.6115 ];
    let mut points = Vec::new();
    for i in 0..num_points {
        points.push(Vec3d::new(ary[3*i], ary[3*i+1], ary[3*i+2]));
    }
    // let mut rng = rand::thread_rng();
    // let points = (0..num_points)
    //     .into_iter()
    //     .map(|_| {
    //         let x = rng.gen();
    //         let y = rng.gen();
    //         let z = x * x + y * y;

    //         Vec3d::new(x, y, z)
    //     })
    //     .collect::<Vec<Vec3d>>();
    let points = sort_and_deduplicate(points);
    let r3s = from_vec3d_to_r3(points);

    println!("point = {:?}", r3s);


    let ts = del3d_rs::newton_apple_delauney(&r3s);

    println!("point = {:?}", ts);
}
