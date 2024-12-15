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
const R_E: f64 = 2.8179403262e-13;  // Classical electron radius in cm

lazy_static! {
    static ref G_A_GAMMA: f64 = 1e-7;         // Axion-photon coupling constant
    static ref TORSION_SCALAR: f64 = 1e-4;   // Torsion strength (arbitrary scaling)
    static ref PHOTON_TO_NEUTRINO_COEFF: f64 = 1e-10; // Photon to neutrino conversion efficiency
}

//----------------------------------------------
// SIMULATION PARAMETERS
//----------------------------------------------
const NX: usize = 7;               // Lattice size
const NY: usize = 7;
const NZ: usize = 7;
const DX: f64 = 0.5e-15;            // Spatial resolution (m)
const DT: f64 = 1.22e-17; //22           // Time step (s)
const STEPS: usize = 100;          // Total simulation steps

// Particle densities
const PHOTON_INIT: f64 = 1e38;
const AXION_INIT: f64 = 1e32;
const NEUTRINO_INIT: f64 = 1e35;

//----------------------------------------------
// FIELD STRUCTURE
//----------------------------------------------
struct Field {
    photon_density: Vec<f64>,
    axion_density: Vec<f64>,
    neutrino_density: Vec<f64>,
    torsion: Vec<f64>,
}

impl Field {
    fn new() -> Self {
        let size = NX * NY * NZ;
        Field {
            photon_density: vec![PHOTON_INIT; size],
            axion_density: vec![AXION_INIT; size],
            neutrino_density: vec![NEUTRINO_INIT; size],
            torsion: vec![*TORSION_SCALAR; size],
        }
    }

    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + NX * (y + NY * z)
    }
}

//----------------------------------------------
// MAIN TIME EVOLUTION
//----------------------------------------------
fn main() {
    let mut field = Field::new();
    let mut file = File::create("results_with_torsion.csv").unwrap();
    writeln!(file, "time(s),avg_photon_density,avg_axion_density,avg_neutrino_density").unwrap();

    for step in 0..STEPS {
        let t = step as f64 * DT;

        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx = field.idx(x, y, z);

                    // Get densities
                    let n_ph = field.photon_density[idx];
                    let n_ax = field.axion_density[idx];
                    let n_nu = field.neutrino_density[idx];
                    let torsion = field.torsion[idx];

                    // Axion to photon conversion
                    let d_ax_to_ph = *G_A_GAMMA * n_ax * DT;

                    // Photon to neutrino conversion via torsion
                    let d_ph_to_nu = *PHOTON_TO_NEUTRINO_COEFF * n_ph * torsion * DT;

                    // Update fields
                    field.photon_density[idx] -= d_ph_to_nu;
                    field.photon_density[idx] += d_ax_to_ph;
                    field.axion_density[idx] -= d_ax_to_ph;
                    field.neutrino_density[idx] += d_ph_to_nu;

                    // Prevent negative densities
                    if field.photon_density[idx] < 0.0 {
                        field.photon_density[idx] = 0.0;
                    }
                    if field.axion_density[idx] < 0.0 {
                        field.axion_density[idx] = 0.0;
                    }
                    if field.neutrino_density[idx] < 0.0 {
                        field.neutrino_density[idx] = 0.0;
                    }
                }
            }
        }

        // Calculate averages
        let vol = (NX * NY * NZ) as f64;
        let avg_photon = field.photon_density.iter().sum::<f64>() / vol;
        let avg_axion = field.axion_density.iter().sum::<f64>() / vol;
        let avg_neutrino = field.neutrino_density.iter().sum::<f64>() / vol;

        writeln!(file, "{},{},{},{}", t, avg_photon, avg_axion, avg_neutrino).unwrap();
    }

    println!("Simulation complete. Results saved to results_with_torsion.csv");
}
