use std::f64::consts::PI;
use std::io::Write;
use lazy_static::lazy_static;

const C: f64 = 2.99792458e8;   // m/s
const H: f64 = 6.62607015e-34; // J·s
const K_B: f64 = 1.380649e-23; // J/K
const G: f64 = 6.67430e-11;    // m^3 kg^-1 s^-2
const TWO: f64 = 2.0;

const G_EARTH: f64 = 9.81; // m/s^2
const R0: f64 = 10.0;      // m

// Grid parameters
const NX: usize = 50;
const NY: usize = 50;
const NZ: usize = 50;
const DX: f64 = 0.1;  // spatial step in meters
const DT: f64 = 1e-2; // time step in seconds 1e-3

lazy_static! {
    static ref LAMBDA_QCD: f64 = TWO.powi(3); // placeholder QCD scale
    static ref T_CMB: f64 = 2.7255; // K
}

// 3D arrays for photon and exotic matter densities
struct Field3D {
    photons: Vec<f64>,
    exotic: Vec<f64>,
}

impl Field3D {
    fn new(nx: usize, ny: usize, nz: usize, photon_init: f64) -> Self {
        let size = nx * ny * nz;
        Field3D {
            photons: vec![photon_init; size],
            exotic: vec![0.0; size],
        }
    }

    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + NX * (y + NY * z)
    }

    // Simple finite difference Laplacian in 3D
    fn laplacian(&self, arr: &Vec<f64>, x: usize, y: usize, z: usize) -> f64 {
        let idx = self.idx(x,y,z);
        let val = arr[idx];
        
        let xm = if x > 0 { arr[self.idx(x-1,y,z)] } else { val };
        let xp = if x < NX-1 { arr[self.idx(x+1,y,z)] } else { val };
        let ym = if y > 0 { arr[self.idx(x,y-1,z)] } else { val };
        let yp = if y < NY-1 { arr[self.idx(x,y+1,z)] } else { val };
        let zm = if z > 0 { arr[self.idx(x,y,z-1)] } else { val };
        let zp = if z < NZ-1 { arr[self.idx(x,y,z+1)] } else { val };

        (xp + xm + yp + ym + zp + zm - 6.0*val)/(DX*DX)
    }
}

// Hypothetical QCD correction
fn qcd_correction(n_photon: f64) -> f64 {
    (n_photon.sqrt()) * *LAMBDA_QCD
}

// Hypothetical exotic matter rate
fn exotic_matter_rate(n_photon: f64) -> f64 {
    let correction = qcd_correction(n_photon);
    1e-35 * correction
}

// Background magnetic field effect (symbolic):
// For instance, let the presence of B-field slightly attenuate photons:
fn magnetic_attenuation_rate(b_field: f64, n_photon: f64) -> f64 {
    // Scale attenuation with B and photon density
    1e-10 * b_field * n_photon
}

fn planck_number_density(T: f64) -> f64 {
    let h = H;
    let kb = K_B;
    let c = C;
    let t = T;

    let nu_min = 1.0e7;
    let nu_max = 1.0e15; // should be 1.0e14
    let steps = 10000;
    let dnu = (nu_max - nu_min)/(steps as f64);

    let mut n_ph = 0.0;
    for i in 0..steps {
        let nu = nu_min + (i as f64)*dnu;
        let x = (h*nu)/(kb*t);
        let bose = 1.0/(f64::exp(x)-1.0);
        let dn = (8.0*PI*nu.powi(2)/(c.powi(3))) * bose * dnu;
        n_ph += dn;
    }
    n_ph
}

fn main() {
    let n_photon_init = planck_number_density(*T_CMB);
    println!("Initial CMB photon number density: {} photons/m^3", n_photon_init);

    // Initialize fields
    let mut field = Field3D::new(NX, NY, NZ, n_photon_init);

    // Assume a uniform background magnetic field along z
    let b_field = 1.0; // 1e20; // Tesla, as a placeholder

    // Open file for output
    let mut file = std::fs::File::create("3d_sim_output.csv").unwrap();
    writeln!(file, "time(s), average_n_photon(m^-3), average_n_exotic(m^-3)").unwrap();

    let total_time = 10e-1 + (1.0+C.powi(4)).sqrt() * 0.00000000000000000125; // 10e-1 * C.powi(20); // run for 0.1 seconds for demo
    let steps = (total_time/DT) as usize;

    // We'll do a simple forward Euler update to show the concept
    for step in 0..steps {
        let mut new_photons = field.photons.clone();
        let mut new_exotic = field.exotic.clone();

        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx = field.idx(x,y,z);
                    let n_ph = field.photons[idx];
                    let n_ex = field.exotic[idx];

                    // Simple diffusion-like photon evolution + attenuation:
                    let lap = field.laplacian(&field.photons, x, y, z);
                    let gamma_B = magnetic_attenuation_rate(b_field, n_ph);

                    // Update photon density: (not physically accurate, just a demonstration)
                    // d(n_ph)/dt = D * lap(n_ph) - gamma_B * n_ph
                    let D = 1e-3; // arbitrary diffusion coefficient
                    let dn_ph = D * lap - gamma_B;
                    new_photons[idx] = n_ph + dn_ph * DT;

                    // Exotic matter formation:
                    // d(n_ex)/dt = exotic_matter_rate(n_ph)
                    let dn_ex = exotic_matter_rate(n_ph);
                    new_exotic[idx] = n_ex + dn_ex * DT;
                }
            }
        }

        field.photons = new_photons;
        field.exotic = new_exotic;

        // Compute averages for output:
        let avg_n_ph: f64 = field.photons.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_n_ex: f64 = field.exotic.iter().sum::<f64>()/(NX*NY*NZ) as f64;

        let time = step as f64 * DT;
        writeln!(file, "{},{},{}", time, avg_n_ph, avg_n_ex).unwrap();
    }

    println!("3D simulation completed. Results in 3d_sim_output.csv");
}
