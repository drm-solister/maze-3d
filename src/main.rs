#![allow(dead_code, unused_variables, unused_mut, unused_imports)]

use std::mem;
use rand::Rng;

mod maze;
mod neighborhood;
mod util;
mod cell_auto;

fn main() {
    // maze_demo();
}

fn maze_demo() {

    let dims: (i32, i32, i32) = (51,1,51);
    // let dims = (11, 3, 11);
    let len = (dims.0*dims.1*dims.2) as usize;

    // place a player and goal
    let start = 0;
    let end = len-1;

    // create the grid
    let mut maze = maze::Maze::new(dims, neighborhood::NeighborhoodMethod::Moore);

    // create the maze
    maze.generate_prims();

    // print the grid with maze
    maze.print();

    // solve the maze
    let x = maze.a_star(start, end).unwrap();

    // add the path taken to the maze
    for i in &x {
        if let maze::Block::Empty = maze.blocks[*i as usize] {
            maze.blocks[*i as usize] = maze::Block::Occupied;
        }
    }

    // print the solved maze
    maze.blocks[end] = maze::Block::Goal;
    maze.print();

}


fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

/*
current problems

8/6/2022
prims now works perfectly! although i think its technically possible for the end point to be blocked in, albiet rare
the more i think about it how the hell can i guarantee any spot is open for a potential goal? other than manually making it open
von neumann still broken prob

*/

