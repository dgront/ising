use std::env;
use clap::Parser;
use ising::Ising2D;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use log::info;

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// number of inner MC cycles
    #[arg(long, default_value_t = 10)]
    inner_cycles: usize,

    /// number of outer MC cycles between frames recorded to PNG files
    #[arg(long, default_value_t = 1)]
    outer_per_frame: usize,

    /// stem of the PNG files to be generated
    #[arg(long, default_value_t = String::from("frame"))]
    png_frames: String,

    /// number of outer MC cycles
    #[arg(long, default_value_t = 100)]
    outer_cycles: usize,

    /// temperature in reduced units
    #[arg(short = 't', long = "temperature", default_value_t = 10.0)]
    temperature: f64,

    /// system size - number of spins along one dimension
    #[arg(short = 'l', long = "length", default_value_t = 16)]
    n_side: usize,

    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn main() {
    let args = Args::parse();
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
        if args.verbose { env::set_var("RUST_LOG", "debug"); }
    }
    env_logger::init();

    info!("spins = {}", args.n_side * args.n_side);
    info!("inner-cycles = {}", args.inner_cycles);
    info!("outer-cycles = {}", args.outer_cycles);
    info!("temperature  = {}", args.temperature);

    let mut model = Ising2D::new(args.n_side);
    let mut rng = StdRng::seed_from_u64(12345);
    model.randomize(&mut rng);
    let filename = format!("{}{:06}.png", args.png_frames, 0);
    model.save_png(&filename).unwrap();

    let mut iframe = 0;
    for outer in 0..args.outer_cycles {
        let mut acceptance = 0.0;
        for _inner in 0..args.inner_cycles {
            acceptance += model.metropolis_sweep_naive(1.0/args.temperature, &mut rng);
        }
        acceptance /= args.inner_cycles as f64;

        println!("{:4} {:.4} {:.4}", outer + 1, model.energy() / (model.n_spins() as f64), acceptance);
        // model.show();
        if (outer + 1) % args.outer_per_frame == 0 {
            iframe += 1;
            let filename = format!("{}{:06}.png", args.png_frames, iframe);
            model.save_png(&filename).unwrap();
        }
    }
}