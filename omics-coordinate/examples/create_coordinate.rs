use omics_coordinate::Coordinate;
use omics_coordinate::system::One;
use omics_coordinate::system::Zero;

pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let coordinate = Coordinate::<Zero>::try_new("seq0", "+", 0)?;
    println!("{:#}", coordinate);

    let coordinate = Coordinate::<One>::try_new("seq0", "+", 1)?;
    println!("{:#}", coordinate);

    Ok(())
}
