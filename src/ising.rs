use std::env;
use std::time::Instant;
use clap::Parser;
use ising::Ising2D;
use rand::{SeedableRng};
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
    #[arg(long)]
    png_frames: Option<String>,

    /// number of outer MC cycles
    #[arg(long, default_value_t = 100)]
    outer_cycles: usize,

    /// Seed for the random number generator
    #[arg(long)]
    seed: Option<usize>,

    /// temperature in reduced units
    #[arg(short = 't', long = "temperature", default_value_t = 10.0)]
    temperature: f64,

    /// system size - number of spins along one dimension (creates a random conformation)
    #[arg(short = 'l', long, required_unless_present = "in_png")]
    n_side: Option<usize>,

    /// system size - number of spins along one dimension
    #[arg(short = 'i', long, required_unless_present = "n_side")]
    in_png: Option<String>,

    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
        if args.verbose { env::set_var("RUST_LOG", "debug"); }
    }
    env_logger::init();

    let mut rng;
    if let Some(seed) = args.seed {
        rng = StdRng::seed_from_u64(seed as u64);
    } else {
        rng = StdRng::from_os_rng();
    }

    let mut model = match (args.n_side, args.in_png) {
        (Some(n), None) => {
            let mut m = Ising2D::new(n);
            m.randomize(&mut rng);
            m
        },
        (None, Some(path)) => Ising2D::from_png(&path)?,
        _ => unreachable!("clap ensures exactly one option"),
    };

    info!("spins = {}", model.n_spins());
    info!("inner-cycles = {}", args.inner_cycles);
    info!("outer-cycles = {}", args.outer_cycles);
    info!("temperature  = {}", args.temperature);



    if let Some(png_stem) = &args.png_frames {
        let filename = format!("{}{:06}.png", png_stem, 0);
        model.save_png(&filename).unwrap();
    }

    let mut iframe = 0;
    let start = Instant::now();
    for outer in 0..args.outer_cycles {
        let mut acceptance = 0.0;
        for _inner in 0..args.inner_cycles {
            // acceptance += model.metropolis_sweep_naive(1.0/args.temperature, &mut rng);
            acceptance += model.metropolis_sweep_checkerboard(1.0/args.temperature, &mut rng);
        }
        acceptance /= args.inner_cycles as f64;
        let largest_per_spin = model.largest_cluster() as f64 / (model.n_spins() as f64);
        println!("{:4} {:.4} {:.4} {:.4} {:.4}",
                 outer + 1, model.energy() / (model.n_spins() as f64),
                 model.magnetization(), largest_per_spin, acceptance);
        // model.show();
        if let Some(png_stem) = &args.png_frames {
            if (outer + 1) % args.outer_per_frame == 0 {
                iframe += 1;
                let filename = format!("{}{:06}.png", png_stem, iframe);
                model.save_png(&filename).unwrap();
            }
        }
    }
    let elapsed = start.elapsed();
    info!("Elapsed: {:.3} s", elapsed.as_secs_f64());

    Ok(())
}