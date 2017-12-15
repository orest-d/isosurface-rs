extern crate isosurface;
extern crate nalgebra;
use isosurface::*;
use nalgebra::{Vector3};
use std::cmp::{Ordering,PartialOrd};

fn partial_max<T:PartialOrd+Copy>(v:&Vec<T>)->Option<T>{
    if v.len()==0{
        None
    }
        else if v.len()==1{
            Some(v[0])
        }
            else{
                let mut acc = v[0];
                for x in v.iter(){
                    match x.partial_cmp(&acc){
                        Some(Ordering::Greater) => {acc = *x}
                        None => {return None}
                        _ => {}
                    }
                }
                Some(acc)
            }
}

#[test]
fn test_basis_grid(){
    let g = Grid::new_box_basis_grid(
        Vector3::new(0.,0.,0.),
        Vector3::new(1.,0.,0.),
        Vector3::new(0.,1.,0.),
        Vector3::new(0.,0.,1.),
        3,3,3);
    assert_eq!(g.points.len(),3*3*3);
    let i:GIndex = 0;
    assert_eq!(g.points[i as usize].x,0.0);
    //println!("edges {}",g.edges.len());
}

#[test]
fn test_relevant_points(){
    let g = Grid::new_box_basis_grid(
        Vector3::new( -1.0,-1.0,-1.0),
        Vector3::new(  0.1, 0.0, 0.0),
        Vector3::new(  0.0, 0.1, 0.0),
        Vector3::new(  0.0, 0.0, 0.1),
        11,11,11);
    let pv = OnTheFlyPointValues{
        points:&g.points,
        function:&|v:Point|{v.norm()-0.001}
    };
    let rp = (&g).relevant_points(&pv);
    let pos = rp.into_iter().position(|x|{x}).unwrap();
    let selected = &g.points[pos];
    println!("P {} {:?} {}",pos,selected,selected.norm());
    assert!(selected.norm()<=0.1);
}

#[test]
fn test_find_intersection() {
    let p = find_intersection(
        Vector3::new(0.,0.,0.),
        Vector3::new(10.,0.,0.),
        &|v: Point| { v.norm() - 1.0 },
        0.001,
        1e-6
    );
    //println!("inter {:?},{}",&p,p.norm());
    assert!((p-Vector3::new(1.,0.,0.)).norm()<0.001)
}
#[test]
fn test_make_mesh() {
    let g = Grid::new_box_basis_grid(
        Vector3::new(-1.,-1.,-1.),
        Vector3::new(0.1,0.,0.),
        Vector3::new(0.,0.1,0.),
        Vector3::new(0.,0.,0.1),
        21,21,21);
    let function = |v:Point|{v.norm()-1.0};
    let pv = OnTheFlyPointValues{
        points:&g.points,
        function:&function
    };

    let mesh = make_mesh(&g,&pv,&function);
    //        println!("mesh points {:?}",mesh.points);
    //        println!("mesh triangles {:?}",mesh.triangles);
    let d = mesh.points.iter().map(|v|{(v.norm()-1.0).abs()}).collect::<Vec<Scalar>>();
    //println!("mesh points max norm diff {}",super::partial_max(&d).unwrap());
    //println!("mesh points norm diff {:?}",mesh.points.iter().map(|v|{(v.norm()-1.0).abs()}).collect::<Vec<f32>>());
    assert!(mesh.points.iter().all(|v|{(v.norm()-1.0).abs()<0.15}));

}

#[test]
fn test_partial_max(){
    assert_eq!(partial_max(&vec![0.1]),Some(0.1));
    assert_eq!(partial_max(&vec![0.1,0.2]),Some(0.2));
    assert_eq!(partial_max(&vec![0.1,0.2,0.0/0.0]),None);
}
