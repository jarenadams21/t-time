use std::f64::consts::PI;
use std::io::Write;

// Physical constants (remains mostly for reference)
const C: f64 = 2.99792458e8;
const G: f64 = 6.67430e-11;
const HBAR: f64 = 1.054571817e-34;
const K_B: f64 = 1.380649e-23;

// Lattice parameters (3D space)
const NX: usize = 20;
const NY: usize = 20;
const NZ: usize = 20;
const DX: f64 = 1.0;     // Spatial step (m)
const DT: f64 = 0.10;    // Time step (s)
const TOTAL_TIME: f64 = 100.0; // Longer simulation time: 100 s
const STEPS: usize = (TOTAL_TIME / DT) as usize;

// Magnetic field scale
const MU_PRIMED: f64 = 4.0 * PI * 10e7;
const B_0: f64 = 0.0; //1.0 * MU_PRIMED;  // 1 Tesla baseline (placeholder)

// QCD-like parameters
const EPSILON_CRIT: f64 = 1e6;
const DELTA: f64 = 1e5; 

// Axion-photon and neutrino parameters (refined to smaller couplings)
const G_A_GAMMA: f64 = 1e-14; // much smaller coupling
const GAMMA_A: f64 = 1e-9;    // reduced axion decay rate

// Neutrino parameters (weaker interactions)
const D_NU: f64 = 1e-4;   // reduced neutrino diffusion
const LAMBDA_NU: f64 = 1e-14; // extremely small energy sink rate

// Metric parameters representing gravitational waves with LIGO-like scale
// Strain amplitude ~ 10^-21 (typical LIGO detection scale)
// Frequencies ~ 100 Hz and 200 Hz
// Wavenumbers correspond to wavelengths of millions of meters.
const H1: f64 = 1e-21;
const H2: f64 = 5e-22;
const OMEGA_1: f64 = 2.0 * PI * 100.0;   // 100 Hz
const OMEGA_2: f64 = 2.0 * PI * 200.0;   // 200 Hz
const K_1: f64 = 2.0 * PI / 3e6;         // wavelength ~3000 km
const K_2: f64 = 2.0 * PI / 1.5e6;       // wavelength ~1500 km

// Diffusion coefficients (keep photons/energy/axions modest)
const D_PH: f64 = 0.1;   
const D_E: f64 = 0.01;   
const D_A: f64 = 0.01;   

struct FluidField {
    energy: Vec<f64>,
    photon_density: Vec<f64>,
    axion_density: Vec<f64>,
    neutrino_density: Vec<f64>,
}

impl FluidField {
    fn new(epsilon_init: f64, photon_init: f64, axion_init: f64, neutrino_init: f64) -> Self {
        let size = NX * NY * NZ;
        FluidField {
            energy: vec![epsilon_init; size],
            photon_density: vec![photon_init; size],
            axion_density: vec![axion_init; size],
            neutrino_density: vec![neutrino_init; size],
        }
    }

    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + NX * (y + NY * z)
    }

    fn laplacian(&self, arr: &Vec<f64>, x: usize, y: usize, z: usize) -> f64 {
        let val = arr[self.idx(x,y,z)];
        let xp = arr[self.idx((x+1)%NX, y, z)];
        let xm = arr[self.idx((x+NX-1)%NX, y, z)];
        let yp = arr[self.idx(x, (y+1)%NY, z)];
        let ym = arr[self.idx(x, (y+NY-1)%NY, z)];
        let zp = arr[self.idx(x, y, (z+1)%NZ)];
        let zm = arr[self.idx(x, y, (z+NZ-1)%NZ)];
        (xp + xm + yp + ym + zp + zm - 6.0*val)/(DX*DX)
    }
}

// Nonlinear QCD-like EoS:
fn qcd_pressure(epsilon: f64) -> f64 {
    if epsilon > EPSILON_CRIT {
        0.33 * epsilon.powf(4.0) // was 1.2
    } else {
        0.33 * epsilon
    }
}

// QGP fraction
fn qgp_fraction(epsilon: f64) -> f64 {
    0.5 * (1.0 + ((epsilon - EPSILON_CRIT) / DELTA).tanh())
}

// Metric factor for gravitational waves
fn metric_factor(t: f64, x: f64) -> f64 {
    let h_combined = H1 * (OMEGA_1*t - K_1*x).cos() + H2 * (OMEGA_2*t - K_2*x).sin();
    1.0 + h_combined
}

fn main() {
    // Initial conditions in a QGP-like regime with large energy density
    let epsilon_init = 1e7;       // J/m^3
    let photon_init = 1e10;       // photons/m^3
    let axion_init = 1e-2;        // axions/m^3
    let neutrino_init = 0.86;     // neutrinos/m^3

    let mut field = FluidField::new(epsilon_init, photon_init, axion_init, neutrino_init);

    let mut file = std::fs::File::create("fluid_lattice.csv").unwrap();
    writeln!(file, "time(s),avg_energy(J/m^3),avg_photon(m^-3),avg_axion(m^-3),avg_neutrino(m^-3),avg_qgp_fraction").unwrap();

    for step in 0..STEPS {
        let t = step as f64 * DT;
        let mut new_energy = field.energy.clone();
        let mut new_photons = field.photon_density.clone();
        let mut new_axions = field.axion_density.clone();
        let mut new_neutrinos = field.neutrino_density.clone();

        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx = field.idx(x,y,z);

                    let e = field.energy[idx];
                    let n_ph = field.photon_density[idx];
                    let n_a = field.axion_density[idx];
                    let n_nu = field.neutrino_density[idx];

                    let p = qcd_pressure(e);
                    let alpha = metric_factor(t, x as f64 * DX);

                    // Laplacians
                    let lap_e = field.laplacian(&field.energy, x,y,z);
                    let lap_ph = field.laplacian(&field.photon_density, x,y,z);
                    let lap_a = field.laplacian(&field.axion_density, x,y,z);
                    let lap_nu = field.laplacian(&field.neutrino_density, x,y,z);

                    // Diffusion updates
                    let de = D_E * lap_e * DT * alpha;
                    let dn_ph = D_PH * lap_ph * DT * alpha;
                    let dn_a = D_A * lap_a * DT * alpha;
                    let dn_nu = D_NU * lap_nu * DT * alpha;

                    // Axion-photon coupling
                    let d_ph_axion = G_A_GAMMA * n_a * B_0.powi(2) * DT;
                    let d_a_loss = -GAMMA_A * n_a * DT;

                    // Neutrino energy sink
                    let d_e_nu = -LAMBDA_NU * n_nu * DT;

                    // QCD-driven sink (mimic expansion)
                    let d_e_qcd = -p * 1e-3 * DT;

                    // Update fields
                    new_photons[idx] = n_ph + dn_ph + d_ph_axion;
                    new_axions[idx] = n_a + dn_a + d_a_loss;
                    new_neutrinos[idx] = n_nu + dn_nu;
                    new_energy[idx] = e + de + d_e_nu + d_e_qcd;

                    // Prevent negatives
                    if new_photons[idx] < 0.0 { new_photons[idx] = 0.0; }
                    if new_axions[idx] < 0.0 { new_axions[idx] = 0.0; }
                    if new_neutrinos[idx] < 0.0 { new_neutrinos[idx] = 0.0; }
                    if new_energy[idx] < 0.0 { new_energy[idx] = 0.0; }
                }
            }
        }

        field.energy = new_energy;
        field.photon_density = new_photons;
        field.axion_density = new_axions;
        field.neutrino_density = new_neutrinos;

        let avg_e = field.energy.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_ph = field.photon_density.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_a = field.axion_density.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_nu = field.neutrino_density.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_qgp = qgp_fraction(avg_e);

        writeln!(file, "{},{},{},{},{},{}", t, avg_e, avg_ph, avg_a, avg_nu, avg_qgp).unwrap();
    }

    println!("Fluid lattice simulation completed. Results in fluid_lattice.csv");
}
