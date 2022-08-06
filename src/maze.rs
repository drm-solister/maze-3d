use rand::Rng;
use rand::prelude::*;
use rand_chacha::ChaCha20Rng;
use std::cmp::Ordering;
// use std::rc::Rc;
// use std::cell::RefCell;
use colored::Colorize;

use crate::neighborhood;

pub struct Maze {
    pub blocks: Vec<Block>,
    pub dims: (i32, i32, i32),
    pub neighbor_offsets: [Option<i32>; 26],
    pub neighborhood_method: neighborhood::NeighborhoodMethod,
}

impl Maze {

    pub fn new(dims: (i32,i32,i32), neighborhood_method: neighborhood::NeighborhoodMethod) -> Maze {

        let mut blocks = vec![Block::Full; dims.0 as usize*dims.1 as usize*dims.2 as usize];
        let neighbor_offsets = neighborhood::get_offsets(dims, &neighborhood_method);

        Maze {
            blocks,
            dims,
            neighbor_offsets,
            neighborhood_method,
        }
    }

    pub fn print(&self) {
        println!("printing grid");
        let mut layer = 0;
        println!("y = {}", layer);

        for (i, block) in self.blocks.iter().enumerate() {
            if i != 0 {
                if i as i32%self.dims.0 == 0 { print!("\n") }
                if i as i32%(self.dims.0 * self.dims.2) == 0 { 
                    print!("\n");
                    layer+=1;
                    println!("y = {}", layer);
                }
            }
            let show_indicies = false;
            if show_indicies {
                match block {
                    Block::Empty    => print!("[{:^3}]", i),
                    Block::Full     => print!("[ x ]"),
                    Block::Occupied => print!("[ o ]"),
                    Block::Goal     => print!("[ ! ]")
                }
            } else {
                match block {
                    Block::Empty    => print!("[ ]"), // full char █
                    Block::Full     => print!("███"),
                    Block::Occupied => print!("{}", "███".green()),
                    Block::Goal     => print!("[!]")
                }
            }
        }
        println!("");
    }

    // get a 3d coordinate from an index of the vec
    pub fn coord_of_index(&self, index: i32) -> Option<(i32, i32, i32)> {
        if index > self.blocks.len() as i32 { return None}

        let dims = self.dims;
        let x = index%dims.0;
        let z = index%(dims.0*dims.2)/dims.0;
        let y = index/(dims.0*dims.2);

        return Some((x as i32, y as i32, z as i32));
    }

    // get the index of the vec corresponding to a 3d coordinate
    pub fn index_of_coord(&self, coord: (i32, i32, i32)) -> Option<usize> {
        // x + z(x_dim) + y(x_dim*z_dim)
        let index = (coord.0 + coord.2*self.dims.0 + coord.1*(self.dims.0*self.dims.2)) as usize;
        if index > self.blocks.len() { return None}
        return Some(index);
    }

    // prevent wrappinga around the sides
    // if the potential neighbor is one less, you cannot be on the left edge -> index % x_dim == 0
    // if the potential neighbor is one more, you cannot be on the right edge -> index+1 & x_dim == 0
    // if the potential neighbor is x_dim more, you cannot be the bottom edge 
    // if the potential neighbor is x_dim less, you cannot be on the top edge
    // get all the indicies of a specific type of neighbor, based on the maze's neighborhood rules 
    pub fn indicies_of_neighbors(&self, target: i32, look_for: Block) -> Vec<i32> {
        let mut vec = Vec::new();
        for i in self.neighbor_offsets {
            if let None = i { continue ;}

            let new_index = target+i.unwrap() as i32;

            if !(new_index >= self.blocks.len() as i32 || new_index < 0) {
                // the neighbor is in bounds
                let dims = self.dims;
                if target % self.dims.0 == 0 && new_index == target-1 { continue; } // avoid wrapping around
                if (target+1) % self.dims.0 == 0 && new_index == target+1 { continue; }

                let z = target%(dims.0*dims.2)/dims.0;

                if z == 0 && new_index == target-self.dims.0 { continue; }
                if z == self.dims.2-1 && new_index == target+self.dims.0 { continue; }


                if look_for == self.blocks[new_index as usize] {
                    vec.push(new_index);
                }
            }
        }
        return vec;
    }

    /*
    Function that implements the a* algorithm, given the start and end index
    Returns a vector that contains the steps to take from start to end in order

    to improve:
    open_list should be a min heap
    */
    pub fn a_star(&self, start: usize, end: usize) -> Option<Vec<usize>> {
        let dims = self.dims;
        let len = (dims.0 * dims.1 * dims.2) as usize;
        let mut traveled: Vec<usize> = vec![0; len];

        let mut open_list: Vec<usize> = Vec::with_capacity(len);
        let mut closed_list: Vec<usize> = Vec::with_capacity(len);
        let mut gscore: Vec<usize> = vec![usize::MAX; len];
        let mut fscore: Vec<usize> = vec![usize::MAX; len];

        gscore[start] = 0;
        fscore[start] = self.h(start, end);
        

        open_list.push(start);

        while open_list.len() > 0 {
            // this will be improved by making openset a min heap
            // this bad
            // current is the node in open_list having the lowest fScore[] value

            // find index of minimum value in fscore. that index is the node that current should be
            let mut index = open_list[0];
            let mut min = usize::MAX;

            for (i,val) in fscore.iter().enumerate() {
                if val < &min {
                    if open_list.contains(&i) {
                        min = *val;
                        index = i;
                    }
                }
            }

            let current = index;

            if current == end {
                // return Some(traveled);
                return Some(Self::construct_path(traveled, start, end));
            }

            // println!("current: {:?}", current);
            
            // this is a stupid and slow way to remove the current value from the open list
            open_list.retain(|&x| x != current);

            // indicies_of_neighbors only finds neighbors of a specific enum. need to find empty and goal blocks
            let mut neighbors = self.indicies_of_neighbors(current as i32, Block::Empty);
            let mut neighboring_goal = self.indicies_of_neighbors(current as i32, Block::Goal);
            neighbors.append(&mut neighboring_goal);

            for n in neighbors {
                let tentative_gscore = gscore[current] + 1;
                if tentative_gscore < gscore[n as usize] {
                    // traveled.push(n as usize); // this could be wrong
                    // traveled.push(current);
                    traveled[n as usize] = current;
                    gscore[n as usize] = tentative_gscore;
                    fscore[n as usize] = tentative_gscore.saturating_add(self.h(n as usize, end));

                    if !open_list.contains(&(n as usize)) {
                        open_list.push(n as usize);
                    }
                }
            }


        }

        return None;
    }

    // helper function that estimates the cost to reach a goal from the index
    pub fn h(&self, index: usize, end: usize) -> usize {
        let index_coord = self.coord_of_index(index as i32).unwrap();
        let end_coord   = self.coord_of_index(end as i32).unwrap();

        let distances = vec![
            (end_coord.0 - index_coord.0).abs(),
            (end_coord.1 - index_coord.1).abs(),
            (end_coord.2 - index_coord.2).abs(),
        ];

        match self.neighborhood_method {
            neighborhood::NeighborhoodMethod::VonNeumann => {
                return *distances.iter().max().unwrap() as usize;
            },
            neighborhood::NeighborhoodMethod::Moore => {
                return distances.iter().sum::<i32>() as usize;
                // return 0;
            },
            _ => panic!("this shouldnt happen"),
        }
    }

    fn construct_path(traveled: Vec<usize>, start: usize, end: usize) -> Vec<usize> {
        let mut path = Vec::new();
        path.push(end);

        let mut step = end;
        loop {
            let next_index = traveled[step];
            path.push(next_index);
            step = next_index;

            if next_index == start { 
                path.reverse();
                return path
            }
        }
    }

}

#[derive(Debug)]
pub struct Color {
    r: i32,
    g: i32,
    b: i32,
}

impl Color {
    pub fn red() -> Self {
        return Color{ r: 255, g: 0, b: 0}
    }

    pub fn black() -> Self {
        return Color{ r: 0, g: 0, b: 0}
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Block {
    Empty,
    Full,
    Occupied,
    Goal,
}

/*
maze generation algorithms
probably dumb to seperate the code like this but lets worry about that later

currently: prim's
later: ?

prims/maybe all of them could work by having potential walls only on every other tile in every direction
*/
impl Maze {

    /*
    not sure if i can call this real prims
    the dimensions have to be odd for this to work.
    */
    pub fn generate_prims(&mut self) {
        // create a sub-grid of nodes

        /*
        Some(false) -> No edge to vertex
        Some(true)  -> Edge to vertex
        None        -> Not a vertex
        */
        let mut prim_nodes: Vec<Option<bool>> = (0..self.blocks.len())
            .map(|x| {
                if self.coord_of_index(x as i32).unwrap().0%2==0 && self.coord_of_index(x as i32).unwrap().2%2==0 {
                    Some(false)
                } else {
                    None
                }
            })
            .collect();

        // println!("-> {:?}", prim_nodes.len());

        let prim_offsets = neighborhood::get_offsets(self.dims, &neighborhood::NeighborhoodMethod::Prims);

        if self.dims.0 % 2 == 0 || self.dims.1 % 2 == 0 || self.dims.2 % 2 == 0 { panic!("grid dimensions must be odd for maze"); }

        let mut rng = rand::thread_rng();
        // let mut rng = ChaCha20Rng::seed_from_u64(2);

        // generate edges by starting at low indicies and only drawing a single edge in the positive x,y,z direction each. so 3 edges. 
        // this makes for no duplicates and every connection

        let mut edges: Vec<Edge> = Vec::new();
        let mut adj_list: Vec<Vec<&Edge>> = vec![Vec::new(); self.blocks.len()];


        // println!("-> {:?}", indicies_of_prim_neighbors(0, &mut self.blocks, self.dims, prim_offsets));
        // // find all the edges that exist
        for (i, node) in prim_nodes
            .iter()
            .enumerate() 
            .filter(|(i, x)| if let Some(_) = x { true } else { false })
            .map(|(i, x)| (i, x.unwrap())) {

            // println!("node: {:?}", i);
            let node_neighbors = indicies_of_prim_neighbors(i as i32, &mut self.blocks, self.dims, prim_offsets);

            for n in node_neighbors {
                let edge = Edge {
                    to: n,
                    from: i,
                    weight: rng.gen(),
                    connected: false,
                };
                edges.push(edge);
            }
        }

        edges.sort_by(|a,b| a.cmp(b));

        // println!("edges: {:?}", edges);

        // // // create adjacency list with the edges
        for (i, edge) in edges.iter().enumerate() {
            // println!("edge: {:?}", edge.borrow());
            adj_list[edge.to].push(&edges[i]);
            adj_list[edge.from].push(&edges[i]);
        }

        

        let starting_node = 0;
        let ending_node = self.blocks.len() - 1;
        let mut p = starting_node;
        self.blocks[ending_node] = Block::Full;


        let start_uncomment_here = 0;

        // start at one node, p
        // look at its edges that are not already connected. 
        // pick the edge with the least weight
        // draw this edge
        // the edge is now connected
        // p is now this edge
        // loop

        // rename p and next_p to point_a, point_b. they are not "next"
        for i in 0..1000 {
            // println!("i: {:?}", i);
            let mut next_p = p;
            let mut mid_coord = 0;
            let mut draw_mid = false;
            for (i, edge) in edges.iter().enumerate() {
                if prim_nodes[edge.to].unwrap() != prim_nodes[edge.from].unwrap() { // only one of the verts connected by the edge is on the tree

                    p = edge.to;
                    next_p = edge.from;

                    mid_coord = self.mid_coord(edge.to, edge.from);
                    draw_mid = true;
                    

                    if next_p == 140 {
                        let stop_here = 0;
                    }

                    edges.remove(i);

                    break;
                }
            }

            self.blocks[p] = Block::Empty;
            self.blocks[next_p] = Block::Empty;
            // let mid_coord = self.mid_coord(p, next_p);
            // draw_mid = true;
            if draw_mid { self.blocks[mid_coord] = Block::Empty; }

            prim_nodes[p] = Some(true);
            prim_nodes[next_p] = Some(true);

            // let x = 0; // debug breakpoint

            // p = next_p;

            // self.print();

        }



    }

    fn mid_coord(&mut self, start: usize, end: usize) -> usize {
        let start_coord = self.coord_of_index(start as i32).unwrap();
        let end_coord = self.coord_of_index(end as i32).unwrap();

        let avg_x = (start_coord.0 + end_coord.0) / 2;
        let avg_y = (start_coord.1 + end_coord.1) / 2;
        let avg_z = (start_coord.2 + end_coord.2) / 2;

        // if (
        //     (start_coord.0 - end_coord.0).abs() > 2 ||
        //     (start_coord.1 - end_coord.1).abs() > 2 ||
        //     (start_coord.2 - end_coord.2).abs() > 2
        // ) {
        //     return 999;
        // }

        // let mut c = 0;
        // if (start_coord.0 - end_coord.0).abs() >= 2 { c+=1; }
        // if (start_coord.1 - end_coord.1).abs() >= 2 { c+=1; }
        // if (start_coord.2 - end_coord.2).abs() >= 2 { c+=1; }

        // if c>1 { return 999; }

        return self.index_of_coord((avg_x, avg_y, avg_z)).unwrap();
    }

}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Edge {
    to: usize,
    from: usize,
    weight: f64,
    connected: bool,
}

impl Edge {
    fn cmp(&self, b: &Self) -> Ordering {
        if self.weight <= b.weight { 
            return Ordering::Less;
        } else {
            return Ordering::Greater;
        }
    }

    fn connect(&mut self) {
        self.connected = true;
    }
}

// how the hell do i create an edge between all adjacent nodes

// duplicate function except itll use a new neighborhood method. this should be consolidated with the first function,
// the first function should take neighbor offsets as an argument
pub fn indicies_of_prim_neighbors(target: i32, nodes: &mut Vec<Block>, dims: (i32, i32, i32), prim_offsets: [Option<i32>; 26]) -> Vec<usize> {
    let mut vec = Vec::new();
    // let neighbor_offsets = neighborhood::get_offsets(dims, &neighborhood::NeighborhoodMethod::Prims);
    for i in prim_offsets {
        if let None = i { continue ;}

        let new_index = target+i.unwrap() as i32;

        if !(new_index >= nodes.len() as i32 || new_index < 0) {
            // the neighbor is in bounds
            if target % dims.0 == 0 && new_index == target-2 { continue; } // avoid wrapping around
            if (target+1) % dims.0 == 0 && new_index == target+2 { continue; }

            let z = target%(dims.0*dims.2)/dims.0;

            if z == 0 && new_index == target-(dims.0*2) { continue; }
            if z == dims.2-1 && new_index == target+(dims.0*2) { continue; }

            vec.push(new_index as usize);

        }
    }
    return vec;
}

// well i messed up because every node has a weight for all its neighbors, which means the weight from node a -> b != weight b -> a
// it might still be fine because its random but it might also not be fine idk

// the wall can trace itself backwards because each real topological edge is represented by two Edge structs, which is bad

// MRE of my problem: https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=5e4e1f1959f96299a6b313fab3bb4889

// more minimal: https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=5e4e1f1959f96299a6b313fab3bb4889


fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

// could give each edge an id that is its index in edges?
// need to rework so that the maze startes all blocked and then the gen algos clear out space
// for some reason some edges are not always being maked as connected when they should be
// its back tracking too, which shouldnt be possible because there should only be one path youll ever go when you come across a node.
// maybe it only gets marked as connected later for some reason

// rework this probably.
/*
idk why it took me so long to realize but i think that the edges dont have to be marked as connected,
rather the nodes have to be marked wheter they are connected to anything or not.

*/

// impl Maze {
//     pub fn invert(&mut self) {
//         for mut block in self.blocks {
//             match block {
//                 Block::Full => block = Block::Empty,
//                 Block::Empty => block = Block::Full,
//                 _ => (),
//             }
//         }
//     }
// }