// element-wise addition of tuples, which i use as coordinates
pub fn add_tuple (a: (i32, i32, i32), b: (i32, i32, i32) ) {

    let mut c: (i32, i32, i32) = (0,0,0);
    for i in 0..2 {
        c.0 = a.0+b.0;
        c.1 = a.1+b.1;
        c.2 = a.2+b.2
    }

    println!("{:?}", c);

}