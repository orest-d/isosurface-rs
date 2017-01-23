extern crate isosurface;
extern crate nalgebra;
use isosurface::{Grid,OnTheFlyPointValues,make_mesh,Point};
use nalgebra::{Vector3,Norm,Cross,Dot};


fn main() {
    let g = Grid::new_box_basis_grid(
        Vector3{ x:-2.0,y:-2.0,z:-2.0},
        Vector3{x:0.05,y:0.0,z:0.0},
        Vector3{x:0.0,y:0.05,z:0.0},
        Vector3{x:0.0,y:0.0,z:0.05},
        81,81,81);
    //let function = |v:Point|{v.x+0.3};
    //let function = |v:Point|{v.norm()-1.0};
    let function = |v:Point|{
        let w=Vector3{ x:v.x+0.3*(5.0*v.y).sin(),y:v.y,z:v.z+0.2*(8.0*v.y).sin()};
        (w.norm()-1.0)*0.5
    };
    let pv = OnTheFlyPointValues{
        points:&g.points,
        function:&function
    };

    let mesh = make_mesh(&g,&pv,&function);
    mesh.write_stl_text("hello.stl");
}
