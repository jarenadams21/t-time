use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::Rng;
use std::f64::consts::PI;

/// Physical constants
const G: f64 = 6.67430e-11;           // Gravitational constant (m^3 kg^-1 s^-2)
const C: f64 = 3.0e8;                 // Speed of light (m/s)
const HBAR: f64 = 1.0545718e-34;      // Reduced Planck constant (J·s)
const K_B: f64 = 1.380649e-23;        // Boltzmann constant (J/K)

// Black hole parameters
const RS: f64 = 1e-6; // Schwarzschild radius (for example)
fn black_hole_mass(rs: f64) -> f64 {
    (rs * C.powi(2)) / (2.0 * G)
}
fn hawking_temperature(m: f64) -> f64 {
    HBAR * C.powi(3) / (8.0 * PI * G * m * K_B)
}

// Lattice parameters
const LATTICE_SIZE: usize = 20;
const N_SWEEPS: usize = 1000;
const DIM: u32 = 3; // 3D lattice
const COUPLING: f64 = 1.0; // coupling constant for field interactions
const MASS_SQ: f64 = 1.0;  // mass^2 term for the scalar field

/// We consider a simple scalar field lattice model:
/// Hamiltonian (discretized) ~ sum over neighbors (phi_x - phi_y)^2 + MASS_SQ * phi_x^2
/// This represents a simple free (or slightly interacting) scalar field.
///
/// We'll use a simple Metropolis algorithm at thermal equilibrium:
/// Probability ~ exp(-H/k_B T)
///
/// Boundary conditions: periodic for simplicity.
/// This does not violate relativity. We are just sampling field configurations at a given temperature.

struct Lattice {
    size: usize,
    field: Vec<f64>,
    temperature: f64,
    rng: StdRng,
}

impl Lattice {
    fn new(size: usize, temperature: f64) -> Self {
        let mut rng = StdRng::seed_from_u64(42);
        let field = (0..size.pow(DIM))
            .map(|_| rng.gen_range(-0.1..0.1)) // small random initial field
            .collect();
        Lattice { size, field, temperature, rng }
    }

    fn index(&self, x: usize, y: usize, z: usize) -> usize {
        x + self.size * (y + self.size * z)
    }

    fn neighbors(&self, x: usize, y: usize, z: usize) -> [(usize,usize,usize); 6] {
        let xm = if x == 0 { self.size-1 } else { x-1 };
        let xp = if x == self.size-1 { 0 } else { x+1 };
        let ym = if y == 0 { self.size-1 } else { y-1 };
        let yp = if y == self.size-1 { 0 } else { y+1 };
        let zm = if z == 0 { self.size-1 } else { z-1 };
        let zp = if z == self.size-1 { 0 } else { z+1 };

        [(xp,y,z),(xm,y,z),(x,yp,z),(x,ym,z),(x,y,zp),(x,y,zm)]
    }

    fn local_energy(&self, x: usize, y: usize, z: usize) -> f64 {
        // Local contribution: (1/2)*sum_neighbors (phi_x - phi_n)^2 + (MASS_SQ/2)*phi_x^2
        let idx = self.index(x,y,z);
        let phi = self.field[idx];
        let mut e = 0.5 * MASS_SQ * phi*phi;

        let neigh = self.neighbors(x,y,z);
        for &(nx,ny,nz) in &neigh {
            let nidx = self.index(nx,ny,nz);
            let d = phi - self.field[nidx];
            e += 0.5 * d*d;
        }
        e
    }

    fn sweep(&mut self) {
        // Metropolis updates
        for _ in 0..self.size.pow(DIM) {
            let x = self.rng.gen_range(0..self.size);
            let y = self.rng.gen_range(0..self.size);
            let z = self.rng.gen_range(0..self.size);
            let idx = self.index(x,y,z);

            let old_phi = self.field[idx];
            let old_e = self.local_energy(x,y,z);

            let new_phi = old_phi + self.rng.gen_range(-0.1..0.1);
            self.field[idx] = new_phi;
            let new_e = self.local_energy(x,y,z);

            let dE = new_e - old_e;
            if dE > 0.0 {
                let prob = (-dE/(K_B * self.temperature)).exp();
                if self.rng.gen::<f64>() > prob {
                    // reject
                    self.field[idx] = old_phi;
                }
            }
        }
    }

    fn measure_energy(&self) -> f64 {
        let mut E = 0.0;
        for x in 0..self.size {
            for y in 0..self.size {
                for z in 0..self.size {
                    // Each local energy counts neighbor pairs twice, but we do not double count if careful:
                    // We'll just sum local_energy and divide by 2 since each bond counted twice.
                    E += self.local_energy(x,y,z);
                }
            }
        }
        E / 2.0
    }

    fn run(&mut self) {
        for sweep in 0..N_SWEEPS {
            self.sweep();
            if sweep % 100 == 0 {
                let e = self.measure_energy();
                println!("Sweep: {}, Energy per site: {}", sweep, e/(self.size.pow(DIM) as f64));
            }
        }
    }
}

fn main() {
    let mass = black_hole_mass(RS);
    let T_hawk = hawking_temperature(mass);

    println!("Black hole mass: {} kg", mass);
    println!("Hawking temperature: {} K", T_hawk);

    // Use Hawking temperature as the system temperature
    let mut lattice = Lattice::new(LATTICE_SIZE, T_hawk);
    lattice.run();

    let final_energy = lattice.measure_energy();
    println!("Final energy per site: {}", final_energy/(LATTICE_SIZE.pow(DIM) as f64));
    println!("Simulation complete with pure statistical mechanics initialization from Hawking radiation temperature. No violations of Special Relativity introduced.");
}
