use image::{GrayImage, Luma};
use rand::Rng;

pub struct Ising2D {
    /// The number of spins along one dimension. The total number of spins is `n_side * n_side`.
    n_side: usize,
    /// The spins of the system, stored in a 1D row-major vector. The spin at position (x, y) can be accessed as `spins[y * n_side + x]`.
    spins: Vec<i8>,
    /// The neighbour above, with periodic boundary conditions. For example, for a 3x3 system, `prev` would be `[1, 2, 0]`.
    prev: Vec<usize>,
    /// The neighbour on the right, with periodic boundary conditions. For example, for a 3x3 system, `next` would be `[2, 0, 1]`.
    next: Vec<usize>,
}

impl Ising2D {
    pub fn new(n_side: usize) -> Self {
        let spins = vec![1; n_side * n_side];
        let prev = (0..n_side).map(|i| (i + n_side - 1) % n_side).collect();
        let next = (0..n_side).map(|i| (i + 1) % n_side).collect();

        Ising2D { n_side, spins, prev, next  }
    }

    /// Loads a 2D Ising model from a PNG file.
    ///
    /// The PNG file should be a grayscale image where white and black pixels represent spin +1 and -1, respectively.
    /// The size of the model is determined by the dimensions of the image, which must be square.
    pub fn from_png(filename: &str) -> Result<Self, image::ImageError> {

        let img = image::open(filename)?.to_luma8();

        assert_eq!(img.width(), img.height(), "Expected a square lattice");

        let n_side = img.width() as usize;
        let mut model = Self::new(n_side);
        model.spins = img.as_raw().iter().map(|&v| if v > 127 { 1 } else { -1 }).collect();

        Ok(model)
    }

    pub fn randomize<R: Rng>(&mut self, rng: &mut R) {
        for s in self.spins.iter_mut() {
            *s = if rng.random_bool(0.5) { 1 } else { -1 };
        }
    }

    pub fn n_spins(&self) -> usize { self.n_side *self.n_side }

    pub fn show(&self) {
        for y in 0..self.n_side {
            for x in 0..self.n_side {
                let s = self.spins[y * self.n_side + x];

                let c = if s > 0 { '#' } else { ' ' };
                print!("{c}");
            }
            println!();
        }
    }

    pub fn save_png(&self, filename: &str) -> Result<(), image::ImageError> {
        let mut img = GrayImage::new(self.n_side as u32, self.n_side as u32);

        for y in 0..self.n_side {
            for x in 0..self.n_side {
                let s = self.spins[y * self.n_side + x];
                let value = if s > 0 { 255 } else { 0 };
                img.put_pixel(x as u32, y as u32, Luma([value]), );
            }
        }

        img.save(filename)
    }

    /// Computes the per-spin magnetization of the Ising system
    pub fn magnetization(&self) -> f64 {
        let m: i64 = self.spins.iter().map(|&s| s as i64).sum();
        return m as f64 / self.n_spins() as f64;
    }

    /// Largest connected cluster of equal spins.
    ///
    /// The method uses an iterative DFS with an explicit stack;
    /// assumes 4-neighbor connectivity and periodic boundaries.
    pub fn largest_cluster(&self) -> usize {

        let n = self.spins.len();
        let mut visited = vec![false; n];
        let mut stack = Vec::new();
        let mut largest = 0usize;

        for start in 0..n {
            if visited[start] { continue; }

            let spin = self.spins[start];
            visited[start] = true;
            stack.push(start);

            let mut cluster_size = 0;

            while let Some(i) = stack.pop() {
                cluster_size += 1;

                let x = i % self.n_side;
                let y = i / self.n_side;

                let neighbors = [
                    y * self.n_side + self.next[x], // right
                    y * self.n_side + self.prev[x], // left
                    self.next[y] * self.n_side + x, // down
                    self.prev[y] * self.n_side + x, // up
                ];

                for j in neighbors {
                    if !visited[j] && self.spins[j] == spin {
                        visited[j] = true;
                        stack.push(j);
                    }
                }
            }
            largest = largest.max(cluster_size);
        }
        return largest;
    }


    /// Computes the total energy of the Ising system
    pub fn energy(&self) -> f64 {
        let mut e: i64 = 0;

        for y in 0..self.n_side {
            let yp = if y + 1 == self.n_side { 0 } else { y + 1 };

            for x in 0..self.n_side {
                let xp = if x + 1 == self.n_side { 0 } else { x + 1 };

                let i = y * self.n_side + x;
                let s = self.spins[i] as i64;

                let right = self.spins[y * self.n_side + xp] as i64;
                let down  = self.spins[yp * self.n_side + x] as i64;

                e -= s * right;
                e -= s * down;
            }
        }

        return e as f64;
    }

    pub fn metropolis_sweep_naive<R: Rng>(&mut self, beta: f64, rng: &mut R) -> f64{
        let n = self.n_side;
        let mut n_accepted = 0;
        for _ in 0..self.n_spins() {
            let i = rng.random_range(0..self.n_spins());

            let x = i % n;
            let y = i / n;

            let left = if x == 0 { n - 1 } else { x - 1 };
            let right = if x + 1 == n { 0 } else { x + 1 };
            let above = if y == 0 { n - 1 } else { y - 1 };
            let below = if y + 1 == n { 0 } else { y + 1 };

            let s = self.spins[i] as i32;

            let neighbour_sum =
                self.spins[y * n + left] as i32
                    + self.spins[y * n + right] as i32
                    + self.spins[above * n + x] as i32
                    + self.spins[below * n + x] as i32;

            let delta_e = 2 * s * neighbour_sum;
            if delta_e <= 0 || rng.random::<f64>() < (-beta * delta_e as f64).exp() {
                self.spins[i] = -self.spins[i];
                n_accepted += 1;
            }
        }
        n_accepted as f64 / self.n_spins() as f64
    }

    pub fn metropolis_sweep_checkerboard<R: Rng>(&mut self, beta: f64, rng: &mut R) -> f64 {

        let mut n_accepted = 0;
        let p4 = (-4.0 * beta).exp();
        let p8 = (-8.0 * beta).exp();
        for parity in 0..2 {
            for y in 0..self.n_side {
                let x0 = (y + parity) & 1;
                let row  = y * self.n_side;
                let up   = self.prev[y] * self.n_side;
                let down = self.next[y] * self.n_side;

                for x in (x0..self.n_side).step_by(2) {
                    let i = row + x;

                    let s = self.spins[i];

                    let nn =
                        self.spins[row + self.prev[x]] // left
                            + self.spins[row + self.next[x]] // right
                            + self.spins[up + x] // above
                            + self.spins[down + x]; // below

                    let delta_e = 2 * s as i32 * nn as i32;

                    let accept = match delta_e {
                        d if d <= 0 => true,
                        4 => rng.random::<f64>() < p4,
                        8 => rng.random::<f64>() < p8,
                        _ => unreachable!(),
                    };

                    if accept {
                        self.spins[i] = -s;
                        n_accepted += 1;
                    }
                }
            }
        }
        n_accepted as f64 / self.n_spins() as f64
    }
}
