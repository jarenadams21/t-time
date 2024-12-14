use std::f64::consts::PI;
use std::io::Write;

// Physical constants
const C: f64 = 2.99792458e8;     // m/s
const H: f64 = 6.62607015e-34;   // J·s
const K_B: f64 = 1.380649e-23;   // J/K
const G: f64 = 6.67430e-11;      // m^3·kg^-1·s^-2
const ALPHA: f64 = 1.0/137.035999084; // Fine structure constant ~1/137
const SIGMA_T: f64 = 6.6524587321e-29; // Thomson cross section (m^2)
const M_E: f64 = 9.10938356e-31; // electron mass (kg)
const E_CHARGE: f64 = 1.602176634e-19; // elementary charge (C)

// Cosmological parameters for recombination era (approx):
// Scale factor at recombination (z ~ 1100):
const A_REC: f64 = 1.0/1100.0;
// CMB temperature at recombination ~3000 K
const T_REC: f64 = 3000.0; 

// Axion-photon coupling upper limit (rough):
// g_{aγ} < 10^-11 GeV^-1 ~ 10^-20 J^-1 for demonstration
// Just an extremely small number to show negligible effect:
const G_AGAMMA: f64 = 1e-20; // J^-1 (approximate order)

// Cosmic magnetic field upper limit ~1e-9 T:
const B_FIELD: f64 = 1e-9;

// Electron density at recombination (approx.):
// Just after recombination: n_e might be around 10^6 m^-3
// This is a rough estimate:
const N_E: f64 = 1e6; 

// Grid parameters (small for demonstration):
const NX: usize = 50;
const NY: usize = 50;
const NZ: usize = 50;
const DX: f64 = 1.0;    // 1 meter cells for demonstration (not realistic)
const DT: f64 = 1e-9;   // small timestep in seconds

struct Field3D {
    photons: Vec<f64>,
    exotic: Vec<f64>,  // Axion-like particle number density
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

    // 3D Laplacian with periodic boundary conditions
    fn laplacian(&self, arr: &Vec<f64>, x: usize, y: usize, z: usize) -> f64 {
        let val = arr[self.idx(x,y,z)];

        let xm = arr[self.idx((x+NX-1)%NX,y,z)];
        let xp = arr[self.idx((x+1)%NX,y,z)];
        let ym = arr[self.idx(x,(y+NY-1)%NY,z)];
        let yp = arr[self.idx(x,(y+1)%NY,z)];
        let zm = arr[self.idx(x,y,(z+NZ-1)%NZ)];
        let zp = arr[self.idx(x,y,(z+1)%NZ)];

        (xp + xm + yp + ym + zp + zm - 6.0*val)/(DX*DX)
    }
}

// Planck distribution approximation for photon number density at temperature T
fn planck_number_density(T: f64) -> f64 {
    // Integrate number density of photons from Planck's law:
    // n_ph = ∫ (8πν²/c³) [1/(exp(hν/kT)-1)] dν from 0 to ∞
    // Known result: n_ph = 20.28 * (T[K])^3 / cm^3 for CMB at low frequencies
    // Convert to SI: n_ph ≈ 2.03e8 * (T/K)^3 m^-3
    // This is a known standard result: n_ph = (16π (k_B T)^3) / (c^3 h^3) ζ(3)
    // ζ(3)≈1.2020569. Let's just compute directly:
    let zeta3 = 1.202056903159594;
    let factor = (16.0 * PI * (K_B * T).powi(3) * zeta3) / (C.powi(3) * H.powi(3));
    factor
}

// Axion-photon conversion rate:
// For simplicity, assume a small conversion probability per unit time depending on B, g_{aγ}.
// Realistically, the conversion requires coherence length L. Set L ~ DX for demonstration.
// Rate ~ (g_{aγ} * B)^2 * c / L  (this is order-of-magnitude, from P ~ ((g_{aγ} B L)/c)^2)
// If L=DX, rate ~ (g_{aγ}^2 * B^2 * C / DX). This is extremely small.
fn axion_conversion_rate() -> f64 {
    let rate = (G_AGAMMA.powi(2) * B_FIELD.powi(2) * C) / DX;
    rate
}

// Diffusion coefficient from Thomson scattering:
// D ~ c * l / 3, l=1/(n_e σ_T)
// We'll scale it down for computational tractability:
fn diffusion_coefficient() -> f64 {
    let mean_free_path = 1.0/(N_E*SIGMA_T); 
    let d = C * mean_free_path / 3.0; 
    // This is enormous (~ 1.5e14 m?), to make simulation stable on meter scale, reduce:
    // We pretend our box represents a scaled version. Let’s just use the real D and accept that
    // in a 50 m box this doesn't represent the full Universe. Physically realistic scaling is huge.
    d
}

// Metric factor: For a FLRW metric, g_{μν}=diag(1,-a²,-a²,-a²).
// Photon number density scales as 1/a³. If we fix a(t)=A_REC for short timescale simulation:
fn scale_factor() -> f64 {
    A_REC
}

fn main() {
    let a = scale_factor();
    let photon_init = planck_number_density(T_REC) / a.powi(3);

    println!("Initial CMB photon number density at recombination: {:.3e} photons/m^3", photon_init);

    let mut field = Field3D::new(NX, NY, NZ, photon_init);

    let axion_rate = axion_conversion_rate();
    let D = diffusion_coefficient();

    let mut file = std::fs::File::create("3d_sim_output.csv").unwrap();
    writeln!(file, "time(s), average_n_photon(m^-3), average_n_exotic(m^-3)").unwrap();

    let total_time = 1e-5; // simulate a very short time
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

                    // Photon diffusion:
                    let lap = field.laplacian(&field.photons, x, y, z);
                    let dn_ph = D * lap * DT;
                    let new_n_ph = n_ph + dn_ph;

                    // Axion (exotic matter) formation: extremely small
                    // d(n_ex)/dt ~ n_ph * axion_rate
                    // This is a gross simplification; in reality it's more complicated.
                    // We show it is negligible:
                    let dn_ex = n_ph * axion_rate * DT;
                    let new_n_ex = n_ex + dn_ex;

                    new_photons[idx] = new_n_ph;
                    new_exotic[idx] = new_n_ex;
                }
            }
        }

        field.photons = new_photons;
        field.exotic = new_exotic;

        let avg_n_ph = field.photons.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_n_ex = field.exotic.iter().sum::<f64>()/(NX*NY*NZ) as f64;

        let time = step as f64 * DT;
        writeln!(file, "{},{},{}", time, avg_n_ph, avg_n_ex).unwrap();
    }

    println!("3D simulation completed. Results in 3d_sim_output.csv");
}
