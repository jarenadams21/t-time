use std::fs::File;
use std::io::Write;
use rand::Rng;
use std::f64::consts::PI;
use lazy_static::lazy_static;

//----------------------------------------------
// PHYSICAL CONSTANTS AND UNITS
//----------------------------------------------
const C: f64 = 2.99792458e8;        // Speed of light (m/s)
const HBAR: f64 = 1.054571817e-34;  // Reduced Planck constant (J·s)
const G: f64 = 6.67430e-11;         // Gravitational constant (m^3/kg/s^2)
const K_B: f64 = 1.380649e-23;      // Boltzmann constant (J/K)

//----------------------------------------------
// COUPLING CONSTANTS
//----------------------------------------------
lazy_static! {
    static ref G_A_GAMMA: f64 = 1e-9;          // Reduced axion-photon coupling to avoid huge flux
    static ref PHOTON_TO_NEUTRINO_COEFF: f64 = 1e-11; // Photon->Neutrino conversion rate reduced
}

//----------------------------------------------
// SIMULATION PARAMETERS
//----------------------------------------------
// We choose smaller densities to avoid overflow:
const NX: usize = 20;
const NY: usize = 20;
const NZ: usize = 20;

const DX: f64 = 0.5e-15;   // 0.5 fm
const DT: f64 = 1.67e-25;  // timestep ~ typical QGP evolution scale
const STEPS: usize = 100;

// Initial densities (scaled down)
const PHOTON_INIT: f64 = 1e33;
const AXION_INIT: f64 = 1e29;
const NEUTRINO_INIT: f64 = 1e31;

// Energy scale ~2 GeV/fm³ ~3.2e35 J/m³
const ENERGY_INIT: f64 = 3.2e35;
const EPSILON_CRIT: f64 = 1.6e35; 
const DELTA: f64 = 0.2e35; 

// Diffusion coefficients
const D_PH: f64 = 5e-4;
const D_AX: f64 = 5e-4;
const D_NU: f64 = 5e-4;
const D_E: f64 = 5e-4;

// Expansion and neutrino sink terms
const LAMBDA_NU: f64 = 1e-6;
const ALPHA_EXPANSION: f64 = 1e-6;

// Gravitational wave parameters (very small effect)
const GW_STR: f64 = 1e-21;
const GW_FREQ: f64 = 1e3; 

//----------------------------------------------
// FIELD STRUCTURE
//----------------------------------------------
struct Field {
    photon_density: Vec<f64>,
    axion_density: Vec<f64>,
    neutrino_density: Vec<f64>,
    energy_density: Vec<f64>,
}

impl Field {
    fn new() -> Self {
        let size = NX * NY * NZ;
        Field {
            photon_density: vec![PHOTON_INIT; size],
            axion_density: vec![AXION_INIT; size],
            neutrino_density: vec![NEUTRINO_INIT; size],
            energy_density: vec![ENERGY_INIT; size],
        }
    }

    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + NX * (y + NY * z)
    }

    fn boundary_index(x: isize, max: usize) -> usize {
        let mut xx = x;
        if xx < 0 {
            xx = 0;
        } else if xx >= max as isize {
            xx = (max as isize) - 1;
        }
        xx as usize
    }

    fn laplacian(&self, arr: &Vec<f64>, x: usize, y: usize, z: usize) -> f64 {
        let xm = Self::boundary_index(x as isize - 1, NX);
        let xp = Self::boundary_index(x as isize + 1, NX);
        let ym = Self::boundary_index(y as isize - 1, NY);
        let yp = Self::boundary_index(y as isize + 1, NY);
        let zm = Self::boundary_index(z as isize - 1, NZ);
        let zp = Self::boundary_index(z as isize + 1, NZ);

        let c = arr[self.idx(x,y,z)];
        let dx2 = DX*DX;
        (arr[self.idx(xp,y,z)] + arr[self.idx(xm,y,z)]
         + arr[self.idx(x,yp,z)] + arr[self.idx(x,ym,z)]
         + arr[self.idx(x,y,zp)] + arr[self.idx(x,y,zm)] - 6.0*c) / dx2
    }
}

//----------------------------------------------
// EQUATION OF STATE FUNCTIONS
//----------------------------------------------
fn eos_pressure(eps: f64) -> f64 {
    let w_qgp = 0.5 * (1.0 + ((eps - EPSILON_CRIT)/DELTA).tanh());
    let p_qgp = (1.0/3.0)*eps;
    let p_hg = 0.15*eps;
    w_qgp*p_qgp + (1.0 - w_qgp)*p_hg
}

// Minimal metric factor from GW
fn metric_factor(t: f64, x: f64) -> f64 {
    1.0 + GW_STR*(2.0*PI*GW_FREQ*t).sin()*x
}

// A clamp function to avoid Inf/NaN
fn clamp_val(x: f64) -> f64 {
    if x.is_nan() || x.is_infinite() {
        1e40 // Some large finite cap
    } else if x > 1e40 {
        1e40
    } else if x < 0.0 {
        0.0
    } else {
        x
    }
}

//----------------------------------------------
// MAIN TIME EVOLUTION
//----------------------------------------------
fn main() {
    let mut field = Field::new();
    let mut file = File::create("results.csv").unwrap();
    writeln!(file, "time(s),avg_photon_density,avg_axion_density,avg_neutrino_density,avg_energy_density").unwrap();

    for step in 0..STEPS {
        let t = step as f64 * DT;

        let mut new_ph = field.photon_density.clone();
        let mut new_ax = field.axion_density.clone();
        let mut new_nu = field.neutrino_density.clone();
        let mut new_e  = field.energy_density.clone();

        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx = field.idx(x, y, z);

                    let n_ph = field.photon_density[idx];
                    let n_ax = field.axion_density[idx];
                    let n_nu = field.neutrino_density[idx];
                    let eps  = field.energy_density[idx];

                    let lap_ph = field.laplacian(&field.photon_density, x, y, z);
                    let lap_ax = field.laplacian(&field.axion_density, x, y, z);
                    let lap_nu = field.laplacian(&field.neutrino_density, x, y, z);
                    let lap_e  = field.laplacian(&field.energy_density, x, y, z);

                    let p = eos_pressure(eps);
                    let x_pos = x as f64 * DX;
                    let mf = metric_factor(t, x_pos);

                    // Axion <-> Photon conversion
                    let d_ax_to_ph = (*G_A_GAMMA) * n_ax * DT;
                    let d_ph_to_nu = (*PHOTON_TO_NEUTRINO_COEFF) * n_ph * DT;

                    let d_e_nu = LAMBDA_NU * n_nu * DT;
                    let d_e_exp = p * ALPHA_EXPANSION * DT;

                    let ph_new = n_ph + D_PH*lap_ph*DT + d_ax_to_ph - d_ph_to_nu;
                    let ax_new = n_ax + D_AX*lap_ax*DT - d_ax_to_ph;
                    let nu_new = n_nu + D_NU*lap_nu*DT + d_ph_to_nu;
                    let e_new = eps + D_E*lap_e*DT - d_e_nu - d_e_exp;

                    new_ph[idx] = clamp_val(ph_new * mf);
                    new_ax[idx] = clamp_val(ax_new * mf);
                    new_nu[idx] = clamp_val(nu_new * mf);
                    new_e[idx]  = clamp_val(e_new * mf);
                }
            }
        }

        field.photon_density = new_ph;
        field.axion_density = new_ax;
        field.neutrino_density = new_nu;
        field.energy_density = new_e;

        let vol = (NX * NY * NZ) as f64;
        let avg_photon = field.photon_density.iter().sum::<f64>() / vol;
        let avg_axion = field.axion_density.iter().sum::<f64>() / vol;
        let avg_neutrino = field.neutrino_density.iter().sum::<f64>() / vol;
        let avg_energy = field.energy_density.iter().sum::<f64>() / vol;

        writeln!(file, "{},{},{},{},{}", t, avg_photon, avg_axion, avg_neutrino, avg_energy).unwrap();
    }

    println!("Simulation complete. Results saved to results.csv");
}
