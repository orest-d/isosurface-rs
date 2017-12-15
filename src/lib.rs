extern crate nalgebra;

use nalgebra::{Vector3};
use std::collections::HashMap;
use std::io::prelude::*;
use std::fs::File;
use std::cmp::{min,max};

pub type Scalar      = f64;
pub type Point       = Vector3<Scalar>;
pub type Vector      = Vector3<Scalar>;
pub type GIndex      = i32;
pub type Edge        = (GIndex,GIndex);
pub type Function    = Fn(Point)->Scalar;

pub fn normal_edge(edge:Edge)->Edge{
    let (a,b)=edge;
    assert!(a!=b);
    (min(a,b),max(a,b))
}

pub struct Tetrahedron{
    pub points : [GIndex;4]//,
    //pub edges  : [GIndex;6]
}

pub trait PointValues{
    fn value(&self, i:GIndex)->Scalar;
}

pub struct OnTheFlyPointValues<'a>{
    pub points:&'a Vec<Point>,
    pub function:&'a Function
}

impl<'a> PointValues for OnTheFlyPointValues<'a>{
    fn value(&self, i:GIndex)->Scalar{
        (*self.function)((*self.points)[i as usize])
    }
}

impl<'a> PointValues for Vec<Scalar>{
    fn value(&self, i:GIndex)->Scalar{
        (*self)[i as usize]
    }
}

pub fn find_intersection(point1:Point,point2:Point,f:&Fn(Point)->Scalar,precision:Scalar,tiny:Scalar)->Point{
    let mut a  = point1;
    let mut b  = point2;
    let mut fa = f(a);
    let mut fb = f(b);
    if fa==0.0 {return a}
    if fb==0.0 {return b}
    if fa*fb>=0.0{
        println!("ERROR");
        return a+(b-a)*0.5;
    }

    if fa>fb{
        a=point2;
        b=point1;
        let fc = fa;
        fa=fb;
        fb=fc;
    }

    assert!(fa<0.0);
    assert!(fb>0.0);
    while ((b-a).norm() > precision) || (fa - fb).abs()<tiny {
        let c=a+(b-a)*0.5;
        let fc = f(c);
        if fc==0.0{return c}
        if fc<0.0{
            a=c;
            fa=fc;
        }
        else{
            b=c;
            fb=fc;
        }
    }
    a+(b-a)*0.5
}

pub struct EdgePoints<'a>{
    pub grid:&'a Grid,
    pub function:&'a Fn(Point)->Scalar,
    pub points: Vec<Point>,
    pub edge_to_point_index: HashMap<Edge,GIndex>,
    pub precision:Scalar,
    pub tiny:Scalar
}

impl<'a> EdgePoints<'a>{
    pub fn new(g:&'a Grid,f:&'a Fn(Point)->Scalar)->EdgePoints<'a>{
        EdgePoints{
            grid:g,
            function:f,
            points:Vec::new(),
            edge_to_point_index:HashMap::new(),
            precision:0.001,
            tiny:1e-7
        }
    }

    pub fn point_index(&mut self,edge:Edge)->GIndex{
        let e=normal_edge(edge);
        {
            if let Some(i) = self.edge_to_point_index.get(&e) {return *i;}
        }
        let p = find_intersection(self.grid.points[e.0 as usize], self.grid.points[e.1 as usize], self.function, self.precision, self.tiny);
        let i = self.points.len() as GIndex;
        self.edge_to_point_index.insert(e, i);
        self.points.push(p);
        i
    }
}

pub struct Mesh{
    pub points:Vec<Point>,
    pub normals:Vec<Point>,
    pub triangles: Vec<(GIndex,GIndex,GIndex)>
}

impl Mesh{
    pub fn write_stl_text(&self,filename:&str)->std::io::Result<()>{
        let mut f = try!(File::create(filename));
        write!(f,"solid {}\n",filename)?;
        for &(a,b,c) in &self.triangles{
            write!(f,"facet normal {} {} {}\n",0,0,0)?;
            write!(f,"  vertex {} {} {}\n",self.points[a as usize].x,self.points[a as usize].y,self.points[a as usize].z)?;
            write!(f,"  vertex {} {} {}\n",self.points[b as usize].x,self.points[b as usize].y,self.points[b as usize].z)?;
            write!(f,"  vertex {} {} {}\n",self.points[c as usize].x,self.points[c as usize].y,self.points[c as usize].z)?;
            write!(f,"endfacet\n")?;
        }
        write!(f,"endsolid {}",filename)?;
        Ok(())
    }
}

pub fn make_mesh<PV:PointValues>(grid:&Grid,point_values:&PV,function:&Fn(Point)->Scalar)->Box<Mesh>{
    let mut edge_points = EdgePoints::new(&grid,function);
    let mut mesh = Mesh{points:Vec::new(),normals:Vec::new(),triangles:Vec::new()};

    fn cut_mmmp(mesh:&mut Mesh,ep:&mut EdgePoints, a:GIndex,b:GIndex,c:GIndex,d:GIndex){
        let positive_point = ep.grid.points[d as usize];
        let i = ep.point_index((a,d));
        let j = ep.point_index((b,d));
        let k = ep.point_index((c,d));
        let positive_direction = ep.grid.points[d as usize]-ep.points[i as usize];
        let u = ep.points[j as usize]-ep.points[i as usize];
        let v = ep.points[k as usize]-ep.points[i as usize];
        if u.cross(&v).dot(&positive_point) > 0.0{
            mesh.triangles.push((i,j,k))
        }
            else{
                mesh.triangles.push((k,j,i))
            }
    }

    fn cut_mmpp(mesh:&mut Mesh,ep:&mut EdgePoints, a:GIndex,b:GIndex,c:GIndex,d:GIndex){
        let positive_point = ep.grid.points[d as usize];
        let ac = ep.point_index((a,c));
        let ad = ep.point_index((a,d));
        let bc = ep.point_index((b,c));
        let bd = ep.point_index((b,d));
        {
            let positive = ep.grid.points[c as usize] - ep.grid.points[a as usize];
            let u = ep.points[bc as usize] - ep.points[ac as usize];
            let v = ep.points[ad as usize] - ep.points[ac as usize];
            if u.cross(&v).dot(&positive) > 0.0 {
                mesh.triangles.push((ac, bc, ad))
            } else {
                mesh.triangles.push((ad, bc, ac))
            }
        }
        {
            let positive = ep.grid.points[d as usize] - ep.grid.points[b as usize];
            let u = ep.points[bc as usize] - ep.points[bd as usize];
            let v = ep.points[ad as usize] - ep.points[bd as usize];
            if u.cross(&v).dot(&positive) > 0.0 {
                mesh.triangles.push((bd, bc, ad))
            } else {
                mesh.triangles.push((ad, bc, bd))
            }
        }
    }

    fn cut_mppp(mesh:&mut Mesh,ep:&mut EdgePoints, a:GIndex,b:GIndex,c:GIndex,d:GIndex){
        let positive_point = ep.grid.points[a as usize];
        let i = ep.point_index((b,a));
        let j = ep.point_index((c,a));
        let k = ep.point_index((d,a));
        let positive_direction = ep.points[i as usize]-ep.grid.points[a as usize];
        let u = ep.points[j as usize]-ep.points[i as usize];
        let v = ep.points[k as usize]-ep.points[i as usize];
        if u.cross(&v).dot(&positive_point) > 0.0{
            mesh.triangles.push((i,j,k))
        }
        else{
            mesh.triangles.push((k,j,i))
        }
    }

    for t in &grid.tetrahedrons{
        let negative = t.points.into_iter().filter(
            |&i|{point_values.value(*i)<0.}
        ).map(|i|*i).collect::<Vec<GIndex>>();
        let positive = t.points.into_iter().filter(
            |&i|{point_values.value(*i)>=0.}
        ).map(|i|*i).collect::<Vec<GIndex>>();
        assert_eq!(negative.len()+positive.len(),4);
        if negative.len()==1{
//            println!("mppp {:?} {:?}",negative,positive);
            cut_mppp(&mut mesh,&mut edge_points, negative[0], positive[0], positive[1], positive[2]);
        }
        else if negative.len()==2{
//            println!("mmpp {:?} {:?}",negative,positive);
            cut_mmpp(&mut mesh,&mut edge_points, negative[0], negative[1], positive[0], positive[1]);
        }
        else if negative.len()==3{
//            println!("mmmp {:?} {:?}",negative,positive);
            cut_mmmp(&mut mesh,&mut edge_points, negative[0], negative[1], negative[2], positive[0]);
        }
    }
    mesh.points = edge_points.points;
    Box::new(mesh)
}

pub struct Grid{
    pub points:Vec<Point>,
    //pub edges:Vec<Edge>,
    //pub points_to_edges:HashMap<Edge,GIndex>,
    pub tetrahedrons:Vec<Tetrahedron>
}

impl Grid{
    pub fn new()->Grid{
        Grid{
            points          : Vec::new(),
            //edges           : Vec::new(),
            //points_to_edges : HashMap::new(),
            tetrahedrons    : Vec::new(),
        }
    }

    pub fn new_box_basis_grid(origin:Point,b1:Vector,b2:Vector,b3:Vector,size1:GIndex,size2:GIndex,size3:GIndex) -> Grid{
        let n1 = size1 + (1-size1%2);
        let n2 = size2 + (1-size2%2);
        let n3 = size3 + (1-size3%2);
        let mut g = Grid::new();
        let index = |x:GIndex,y:GIndex,z:GIndex|{z+n3*(y+n2*x)};
        let tetrahedron_box = |g:&mut Grid, i:GIndex,ii:GIndex,j:GIndex,jj:GIndex,k:GIndex,kk:GIndex|{
            g.add_tetrahedron(
                index(i, j, k ),
                index(ii,j, k ),
                index(i, jj,k ),
                index(i, j, kk),
            );
            g.add_tetrahedron(
                index(ii,jj,k ),
                index(ii, j,k ),
                index(i, jj,k ),
                index(ii,jj,kk),
            );
            g.add_tetrahedron(
                index(i, jj,kk),
                index(i, j, kk),
                index(ii,jj,kk),
                index(i, jj,k ),
            );
            g.add_tetrahedron(
                index(ii,j, kk),
                index(ii,jj,kk),
                index(i, j, kk),
                index(ii,j, k ),
            );
            g.add_tetrahedron(
                index(ii,jj,kk),
                index(i, j, kk),
                index(ii,j, k ),
                index(i, jj,k ),
            );
        };
        let mut symmetric_box = |g:&mut Grid, i:GIndex,j:GIndex,k:GIndex|{
            tetrahedron_box(g, i   ,i-1 ,j   ,j-1  ,k   ,k-1 );
            tetrahedron_box(g, i-2 ,i-1 ,j   ,j-1  ,k   ,k-1 );
            tetrahedron_box(g, i   ,i-1 ,j-2 ,j-1  ,k   ,k-1 );
            tetrahedron_box(g, i-2 ,i-1 ,j-2 ,j-1  ,k   ,k-1 );
            tetrahedron_box(g, i   ,i-1 ,j   ,j-1  ,k-2 ,k-1 );
            tetrahedron_box(g, i-2 ,i-1 ,j   ,j-1  ,k-2 ,k-1 );
            tetrahedron_box(g, i   ,i-1 ,j-2 ,j-1  ,k-2 ,k-1 );
            tetrahedron_box(g, i-2 ,i-1 ,j-2 ,j-1  ,k-2 ,k-1 );
        };
        for i in 0..n1{
            for j in 0..n2{
                for k in 0..n3{
                    let point = origin + b1*(i as Scalar) + b2*(j as Scalar) + b3*(k as Scalar);
                    //println!("ADD Point {},{},{}: {:?}",i,j,k,point);
                    g.points.push(point);

                    if i>1 && j>1 && k>1 && i%2==0 && j%2==0 && k%2==0 {
                        symmetric_box(&mut g, i,j,k)
                    }
                }
            }
        }
        g
    }
    fn add_tetrahedron(&mut self,a:GIndex,b:GIndex,c:GIndex,d:GIndex){
        let t = Tetrahedron{
            points:[a,b,c,d]
        };
        self.tetrahedrons.push(t);
    }

    pub fn relevant_points<PV:PointValues>(&self,point_values:&PV)->Vec<bool>{
        let mut relevant:Vec<bool> = vec![false;self.points.len()];
        for t in &self.tetrahedrons{
            let v = t.points.into_iter().map(|i|{point_values.value(*i)>0.}).collect::<Vec<bool>>();
            if (&v).into_iter().any(|x|{*x}) && (&v).into_iter().any(|x|{!*x}){

                for point_index in t.points.into_iter(){
                    relevant[*point_index as usize] = true;
                }
            }
        }
        relevant
    }

}
