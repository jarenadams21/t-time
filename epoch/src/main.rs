use std::io::Write;
use lazy_static::lazy_static;
use std::f64::consts::PI;

const C: f64 = 2.99792458e8;   // m/s
const H: f64 = 6.62607015e-34; // J·s
const K_B: f64 = 1.380649e-23; // J/K
const G: f64 = 6.67430e-11;    // m^3 kg^-1 s^-2
const TWO: f64 = 2.0;

// Grid parameters
const NX: usize = 50;
const NY: usize = 50;
const NZ: usize = 50;
const DX: f64 = 0.1;      // spatial step in meters
const DT: f64 = 1e-4;     // time step in seconds

lazy_static! {
    static ref LAMBDA_QCD: f64 = TWO.powi(3); // placeholder QCD scale
    static ref T_CMB: f64 = 2.7255; // K
}

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
fn magnetic_attenuation_rate(b_field: f64, n_photon: f64) -> f64 {
    // Scale attenuation with B and photon density
    1e-10 * b_field * n_photon
}

// Planck distribution approximation for CMB photon number density
fn planck_number_density(T: f64) -> f64 {
    let h = H;
    let kb = K_B;
    let c = C;
    let t = T;

    let nu_min = 1.0e7;
    let nu_max = 1.0e15;
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

// Introduce scaling relevant to black hole holography scenarios:
// Let's consider a black hole mass scale for dimensionless parameters
// (For a real scenario, choose a mass of interest. Here, just an example.)
fn characteristic_scales() -> (f64, f64, f64) {
    let mass_bh = 1.0e30; // kg, roughly solar mass
    let r_s = 2.0 * G * mass_bh / (C*C); // Schwarzschild radius
    let char_time = r_s / C; // characteristic time scale
    (r_s, char_time, mass_bh)
}

// Placeholder metric factor function (could be extended):
fn metric_factor(_x: usize, _y: usize, _z: usize) -> f64 {
    // In a real scenario, this would depend on g_{μν} and possibly coordinates.
    1.0
}

fn main() {
    let n_photon_init = planck_number_density(*T_CMB);
    println!("Initial CMB photon number density: {} photons/m^3", n_photon_init);

    // Initialize fields
    let mut field = Field3D::new(NX, NY, NZ, n_photon_init);

    let b_field = 1.0; // Tesla as placeholder

    // Compute characteristic scales:
    let (r_s, char_time, _mass_bh) = characteristic_scales();

    // Dimensionless diffusion coefficient:
    let dimensionless_diffusion = 1e-24; 
    // Physical D in m^2/s (scaling with black hole radius and time):
    let D = dimensionless_diffusion * (r_s * r_s / char_time);

    // Open file for output
    let mut file = std::fs::File::create("3d_sim_output.csv").unwrap();
    writeln!(file, "time(s), average_n_photon(m^-3), average_n_exotic(m^-3)").unwrap();

    let total_time = 0.1 + (1.0 + C.powi(4)).sqrt() * 0.00000000000000000125;
    let steps = (total_time/DT) as usize;

    for step in 0..steps {
        let mut new_photons = field.photons.clone();
        let mut new_exotic = field.exotic.clone();

        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx = field.idx(x,y,z);
                    let n_ph = field.photons[idx];
                    let n_ex = field.exotic[idx];

                    let alpha_g = metric_factor(x,y,z);

                    let lap = field.laplacian(&field.photons, x, y, z);
                    let gamma_B = magnetic_attenuation_rate(b_field, n_ph) * alpha_g;

                    // Photon evolution incorporating geometric factor:
                    let dn_ph = (D * lap - gamma_B) * alpha_g;
                    let new_n_ph = n_ph + dn_ph * DT;

                    // Exotic matter formation rate scaled by geometry
                    let dn_ex = exotic_matter_rate(n_ph) * alpha_g;
                    let new_n_ex = n_ex + dn_ex * DT;

                    new_photons[idx] = new_n_ph;
                    new_exotic[idx] = new_n_ex;
                }
            }
        }

        field.photons = new_photons;
        field.exotic = new_exotic;

        // Compute averages for output:
        let avg_n_ph = field.photons.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_n_ex = field.exotic.iter().sum::<f64>()/(NX*NY*NZ) as f64;

        let time = step as f64 * DT;
        writeln!(file, "{},{},{}", time, avg_n_ph, avg_n_ex).unwrap();
    }

    println!("3D simulation completed. Results in 3d_sim_output.csv");
}
