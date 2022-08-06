#![allow(dead_code, unused_variables, unused_mut, unused_imports)]

use std::mem;
use rand::Rng;

mod maze;
mod neighborhood;
mod util;

fn main() {

    // let empty_block = maze::Block::Empty(maze::Color::red());
    // let full_block = maze::Block::Full(maze::Color::black());

    let dims: (i32, i32, i32) = (51,1,51);
    // let dims = (11, 3, 11);
    let len = (dims.0*dims.1*dims.2) as usize;
    

    // hard code a maze
    // let wall_indicies = [5,15,35,45];
    // for i in wall_indicies {
    //     blocks[i] = maze::Block::Full;
    // }

    // place a player and goal
    let start = 0;
    let end = len-1;

    // create the grid
    let mut maze = maze::Maze::new(dims, neighborhood::NeighborhoodMethod::Moore);

    // maze.blocks[start] = maze::Block::Occupied;
    // maze.blocks[end] = maze::Block::Goal;

    // maze.print();

    // create the maze
    // maze.generate_maze_hunt_and_kill();
    maze.generate_prims();

    // invert grid cuz im dumb
    // maze.invert();

    // print the grid
    maze.print();

    // solve the maze
    let x = maze.a_star(start, end).unwrap();

    // add the path taken to the maze
    for i in &x {
        if let maze::Block::Empty = maze.blocks[*i as usize] {
            maze.blocks[*i as usize] = maze::Block::Occupied;
        }
    }

    maze.blocks[end] = maze::Block::Goal;
    // print the solved maze
    maze.print();
    // println!("{:?}", x);
}


fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

/*
current problems

the type of coordinates and indicies are different but i think ill have to use both of them in similar contexts
can i make functions that can take either as an argument?
neighbor offsets is an array of length 26 that will have many Nones for Moore. Might want to make this a vec
heuristic function may be slow
reimplement openset as a priority queue or min-heap (*shudders*)


a* kinda going off the shits but unfortunately i do not know why
i calculated the neighbors wrong, it can wrap around the sides oops
it now gets at least one path correct but its going in another branch in the wrong direction. my heuristic may be wrong
just kidding i had to use the construct path function now it seems like it works :flushed:

more testing tomorrow, and in 3d :O

in 3d the neighbors are messed up on bottom and top edges, as opposed to side edges

a* now works 100% i think
with vonneuman rules it goes in a path that sortof looks suboptimal, but it isnt so idk how to fix that really
with vonneuman rules it has new ways to wrap around the edge of the map uuugugghghhh

reference prim's algorithm for maze generation, as well as this github repo: https://github.com/conorpp/3d-maze-generator

i am terrible at writing this maze algorithm
also von neumman rules still wrap the edges

hunt and kill is only going in a straight line now wth

just realized hunt and kill doesnt let me define the start and end point


my incorrect implementation of prims allows for multiple paths, 

8/6/2022
prims now works perfectly! although i think its technically possible for the end point to be blocked in, albiet rare

*/

