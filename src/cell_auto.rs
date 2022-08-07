// cellular automata
// quite a bit of copied and pasted code because its way simpler than trying to seperate prior code into chunks.
// main reason it would be difficult is cell_auto and maze use different enums for the "blocks"/"cells"
// a lot of their functionality does not overlap

use crate::neighborhood;

pub struct CellAuto {
    pub cells: Vec<Cell>,
    pub dims: (i32, i32, i32),
    pub neighbor_offsets: [Option<i32>; 26],
    pub neighborhood_method: neighborhood::NeighborhoodMethod,
}

#[derive(Debug, Clone, Copy)]
pub enum Cell {
    Alive(i16),
    Dead,
}

impl CellAuto {

    pub fn new(dims: (i32,i32,i32), neighborhood_method: neighborhood::NeighborhoodMethod) -> CellAuto {

        let mut cells = vec![Cell::Dead; dims.0 as usize*dims.1 as usize*dims.2 as usize];
        let neighbor_offsets = neighborhood::get_offsets(dims, &neighborhood_method);

        CellAuto {
            cells,
            dims,
            neighbor_offsets,
            neighborhood_method,
        }
    }

    pub fn step(&self) {
        // update every cell based on its neighbors
        // use rayon here
        todo!();
    }
}