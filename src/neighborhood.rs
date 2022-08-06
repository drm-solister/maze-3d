use crate::maze;
use crate::util;

#[derive(Debug)]
pub enum NeighborhoodMethod {
    VonNeumann,
    Moore,
    Diagonal,
    Prims,
}

// all the index offsets that would get you every neighbor
pub fn get_offsets(dims: (i32, i32, i32), neighborhood: &NeighborhoodMethod) -> [Option<i32>; 26] {
    let mut offsets = [None; 26];

    let adj_coords = match neighborhood {
        NeighborhoodMethod::VonNeumann => &VON_NEUMANN_NEIGHBORS[..],
        NeighborhoodMethod::Moore => &MOORE_NEIGHBORS[..],
        NeighborhoodMethod::Diagonal => &DIAGONAL_NEIGHBORS[..],
        NeighborhoodMethod::Prims => &PRIMS_NEIGHBORS[..],
    };

    // println!("====");
    // println!("{:?}", adj_coords);

    for (i,coord) in adj_coords.iter().enumerate() {
        let offset = coord.0 + coord.2*dims.0 + coord.1*(dims.0*dims.2);
        offsets[i] = Some(offset);
    }

    return offsets;
}

pub static VON_NEUMANN_NEIGHBORS: [(i32, i32, i32); 26] = [
    (-1, -1, -1),
    (-1, -1, 0),
    (-1, -1, 1),
    (-1, 0, -1),
    (-1, 0, 0),
    (-1, 0, 1),
    (-1, 1, -1),
    (-1, 1, 0),
    (-1, 1, 1),
    (0, -1, -1),
    (0, -1, 0),
    (0, -1, 1),
    (0, 0, -1),
    (0, 0, 1),
    (0, 1, -1),
    (0, 1, 0),
    (0, 1, 1),
    (1, -1, -1),
    (1, -1, 0),
    (1, -1, 1),
    (1, 0, -1),
    (1, 0, 0),
    (1, 0, 1),
    (1, 1, -1),
    (1, 1, 0),
    (1, 1, 1)
];

pub static MOORE_NEIGHBORS: [(i32, i32, i32); 6] = [
    (-1, 0, 0),
    (1, 0, 0),
    (0, -1, 0),
    (0, 1, 0),
    (0, 0, -1),
    (0, 0, 1)
];

pub static DIAGONAL_NEIGHBORS: [(i32, i32, i32); 20] = [
    (-1, -1, -1),
    (-1, -1, 0),
    (-1, -1, 1),
    (-1, 0, -1),
    // (-1, 0, 0),
    (-1, 0, 1),
    (-1, 1, -1),
    (-1, 1, 0),
    (-1, 1, 1),
    (0, -1, -1),
    // (0, -1, 0),
    (0, -1, 1),
    // (0, 0, -1),
    // (0, 0, 1),
    (0, 1, -1),
    // (0, 1, 0),
    (0, 1, 1),
    (1, -1, -1),
    (1, -1, 0),
    (1, -1, 1),
    (1, 0, -1),
    // (1, 0, 0),
    (1, 0, 1),
    (1, 1, -1),
    (1, 1, 0),
    (1, 1, 1)
];

pub static PRIMS_NEIGHBORS: [(i32, i32, i32); 3] = [
    // (-2, 0, 0),
    (2, 0, 0),
    // (0, -2, 0),
    (0, 2, 0),
    // (0, 0, -2),
    (0, 0, 2)
];