use std::fs::File;
use std::io::Write;
use rand::Rng;
use std::f64::consts::PI;
use lazy_static::lazy_static;

//----------------------------------------------
// PHYSICAL CONSTANTS AND UNITS
//----------------------------------------------
// Fundamental constants (SI Units)
const C: f64 = 2.99792458e8;        // speed of light (m/s)
const HBAR: f64 = 1.054571817e-34;  // reduced Planck (J·s)
const K_B: f64 = 1.380649e-23;      // Boltzmann (J/K)
const G: f64 = 6.67430e-11;         // gravitational constant (m^3/kg/s^2)
const R_E: f64 = 2.8179403262e-13; // classical electron radius in cm
lazy_static! {
    static ref BETA: f64 = 1.0;
    // Thomson cross section:
    static ref SIGMA_T: f64 = (8.0 * PI / 3.0) * R_E.powi(2);
    // Axion: hypothetical low density ~1e32 m^-3
    static ref AXION_INIT: f64 = AXION_INIT_MASS_VAL.powf(1.0);  // The singularity
}
// Conversion factors
// 1 fm = 1e-15 m
// Energy scale: 1 GeV = 1.6021765e-10 J
// Typical QGP energy densities ~1-10 GeV/fm^3 => Convert to J/m^3
// 1 fm^3 = (1e-15 m)^3 = 1e-45 m^3, so 1 GeV/fm^3 = 1.6021765e-10 J / 1e-45 m^3 = 1.6021765e35 J/m^3.
// use a scale: epsilon_init ~ 2 GeV/fm^3 ~ 3.2e35 J/m^3 as starting point.


//----------------------------------------------
// SIMULATION PARAMETERS
//----------------------------------------------
// Lattice size: small box ~ (10 fm)^3 for demonstration
const NX: usize = 10;
const NY: usize = 10;
const NZ: usize = 10;

// Spatial resolution
const DX_FM: f64 = 0.5;          // 0.5 fm per cell
const DX: f64 = DX_FM * 1e-15;   // in meters

// Time step ~0.05 fm/c: 1 fm/c ~ 3.3356e-24 s, so 0.05 fm/c ~1.6678e-25 s
//const DT: f64 = 1.67e-25;
const DT: f64 = 1.22e-22;
// Total time ~ 100 steps => 100*DT ~ 1.67e-23 s (a realistic QGP lifetime)
const STEPS: usize = 1000;

//----------------------------------------------
// QGP PARAMETERS
//----------------------------------------------
// Critical temperature Tc ~170 MeV => critical energy density ~0.5-1 GeV/fm^3
// set EPSILON_CRIT to ~1 GeV/fm^3 ~1.6e35 J/m^3
const EPSILON_CRIT: f64 = 1.6e35;  
const DELTA: f64 = 0.2e35;        // transition width

// Initial energy density ~2 GeV/fm^3
const EPSILON_INIT: f64 = 3.2e35; 

// Photon, Axion, Neutrino initial densities
// Photon: high density in QGP ~ 1e38 m^-3 (arbitrary large number for radiation-dominated plasma)
// Axion: hypothetical low density ~1e32 m^-3
// Neutrino: from thermal estimates ~1e35 m^-3
const PHOTON_INIT: f64 = 1e38;
const AXION_INIT_MASS_VAL: f64 = 1e32;
const NEUTRINO_INIT: f64 = 1e35;


//----------------------------------------------
// NEUTRINO MASS AND STOCHASTICITY
//----------------------------------------------
// Base neutrino mass: 0.15 eV = 0.15 * 1.6021765e-19 J/eV = 2.40326475e-20 J
// We'll just note it, though the code does not explicitly use neutrino mass in computations other than for reference.
// We allow a random mass between 0.05 to 0.8 eV at initialization step if needed for a coupling calculation.
// For simplicity in this code, neutrino mass doesn't directly enter the fluid equations, but could affect neutrino interactions.
const NEUTRINO_MASS_MIN: f64 = 0.05;
const NEUTRINO_MASS_MAX: f64 = 0.8;
fn neutrino_mass_eV() -> f64 {
    let mut rng = rand::thread_rng();
    rng.gen_range(NEUTRINO_MASS_MIN..NEUTRINO_MASS_MAX)
}

//----------------------------------------------
// DIFFUSION AND COUPLINGS
//----------------------------------------------
const D_PH: f64 = 0.05;     // photon diffusion (dimensionless)
const D_E: f64 = 0.01;      // energy diffusion (dimensionless)
const D_A: f64 = 0.01;      // axion diffusion
const D_NU: f64 = 0.01;     // neutrino diffusion

// Axion-photon coupling (small but non-zero)
const G_A_GAMMA: f64 = 1e-7;  
// Axion decay rate
const GAMMA_A: f64 = 1e-3;     

// Neutrino energy sink factor
const LAMBDA_NU: f64 = 1e-5;  

// Introduce a non-zero background magnetic field to enable axion-photon conversions
// For a neutron star interior, fields can be enormous: ~10^11 T
const B_0: f64 = 1e11; 

//----------------------------------------------
// GRAVITATIONAL WAVE METRIC FACTOR
//----------------------------------------------
// Assume frequencies and wavevectors consistent with neutron star mergers
// f ~ kHz range => Let's pick f1=2kHz, f2=4kHz for demonstration
// Strain amplitude ~1e-21
const H1: f64 = 1e-21;
const H2: f64 = 0.5e-21;
const OMEGA_1: f64 = 2.0 * PI * 2000.0;  // rad/s
const OMEGA_2: f64 = 2.0 * PI * 4000.0;  // rad/s
// Wavelength ~1000 m scale is huge compared to fm scale, but we scale down artificially for demonstration
// This is a stretch: we pick k_1, k_2 extremely large to match small domain ~ fm
// Let’s pick k_1 = 1e15 m^-1, k_2=2e15 m^-1 for a short-wavelength GW in nuclear matter
const K_1: f64 = 1e15;
const K_2: f64 = 2e15;

fn metric_factor(t: f64, x: f64) -> f64 {
    let h_combined = H1 * (OMEGA_1*t - K_1*x).cos() + H2 * (OMEGA_2*t - K_2*x).sin();
    1.0 + h_combined
}

//----------------------------------------------
// EQUATION OF STATE
//----------------------------------------------
// We use a piecewise EoS from lattice QCD:
// For epsilon > EPSILON_CRIT, QGP phase: p ~ 1/3 * epsilon (ultrarelativistic limit).
// For epsilon < EPSILON_CRIT, hadron gas: p ~ 0.15 * epsilon.
// Smooth crossover can be mimicked by a tanh blend.

fn qcd_pressure(epsilon: f64) -> f64 {
    let w_qgp = 0.5 * (1.0 + ((epsilon - EPSILON_CRIT) / DELTA).tanh()); 
    let p_qgp = (1.0/3.0)*epsilon;
    let p_hg = 0.15 * epsilon;
    w_qgp * p_qgp + (1.0 - w_qgp)*p_hg
}

fn qgp_fraction(epsilon: f64) -> f64 {
    0.5 * (1.0 + ((epsilon - EPSILON_CRIT) / DELTA).tanh())
}

//----------------------------------------------
// FIELD STRUCTURE
//----------------------------------------------
struct FluidField {
    energy: Vec<f64>,
    photon_density: Vec<f64>,
    axion_density: Vec<f64>,
    neutrino_density: Vec<f64>,
}

impl FluidField {
    fn new() -> Self {
        let size = NX * NY * NZ;
        FluidField {
            energy: vec![EPSILON_INIT; size],
            photon_density: vec![PHOTON_INIT; size],
            axion_density: vec![*AXION_INIT; size],
            neutrino_density: vec![NEUTRINO_INIT; size],
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

//----------------------------------------------
// MAIN TIME EVOLUTION
//----------------------------------------------
fn main() {
    // Random neutrino mass assignment just for record
    let nu_mass = neutrino_mass_eV();

    let mut field = FluidField::new();
    let mut file = File::create("results.csv").unwrap();
    writeln!(file, "time(s),avg_energy(J/m^3),avg_photon(m^-3),avg_axion(m^-3),avg_neutrino(m^-3),avg_qgp_fraction,neutrino_mass_eV").unwrap();

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
                    let dph = D_PH * lap_ph * DT * alpha;
                    let da = D_A * lap_a * DT * alpha;
                    let dnu = D_NU * lap_nu * DT * alpha;

                    // Axion-photon coupling: n_ph increases due to axions in a strong B field
                    let d_ph_axion = G_A_GAMMA * n_a * B_0 * B_0 * DT;
                    let d_a_loss = -GAMMA_A * n_a * DT;

                    // Neutrino energy sink
                    let d_e_nu = -LAMBDA_NU * n_nu * DT;

                    // Expansion-like sink from pressure work (mimic system losing energy)
                    // On fm scale, if we consider a slight expansion dV/V ~ small:
                    // For simplicity:
                    let d_e_qcd = -p * 1e-3 * DT; // small "expansion" factor

                    // Update fields
                    new_photons[idx] = n_ph + dph + d_ph_axion;
                    new_axions[idx] = n_a + da + d_a_loss;
                    new_neutrinos[idx] = n_nu + dnu;
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

        let vol = (NX*NY*NZ) as f64;
        let avg_e = field.energy.iter().sum::<f64>()/vol;
        let avg_ph = field.photon_density.iter().sum::<f64>()/vol;
        let avg_a = field.axion_density.iter().sum::<f64>()/vol;
        let avg_nu = field.neutrino_density.iter().sum::<f64>()/vol;
        let avg_qgp = qgp_fraction(avg_e);

        writeln!(file, "{},{},{},{},{},{},{}", t, avg_e, avg_ph, avg_a, avg_nu, avg_qgp, nu_mass).unwrap();
    }

    println!("Simulation completed. Results in results.csv");
}
