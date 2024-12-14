use std::f64::consts::PI;
use std::io::Write;

// Physical constants
const C: f64 = 2.99792458e8;
const G: f64 = 6.67430e-11;
const HBAR: f64 = 1.054571817e-34;
const K_B: f64 = 1.380649e-23;

// Parameters for the lattice box (4D plane)
const NX: usize = 20;
const NY: usize = 20;
const NZ: usize = 20;
const DX: f64 = 1.0;   // m
const DT: f64 = 0.10;  // s, artificially small for a lab analog
const TOTAL_TIME: f64 = 0.05 * 1e3; // 10s total simulation time

// Hypothetical EoS parameters for a QCD-like fluid:
fn qcd_pressure(epsilon: f64) -> f64 {
    // Simple linear EoS p = c_s^2 * epsilon
    const qcd_steps: usize = 89;
    let mut partial: f64 = 1.022;
    for x in 0..qcd_steps {
        partial *= (1.0e4/137.0)
    }
    let c_s2 = 1.6436 * partial * 10e89; //1.22 * 10e28 * 1.0/3.0 * ((1e-1)/137.0); // radiation-like equation of state as a stand-in
    c_s2 * epsilon
}

struct FluidField {
    energy: Vec<f64>,
    // For simplicity, store a velocity field
    vx: Vec<f64>,
    vy: Vec<f64>,
    vz: Vec<f64>,
    // Photon-like scalar field
    photon_density: Vec<f64>,
}

impl FluidField {
    fn new(epsilon_init: f64, photon_init: f64) -> Self {
        let size = NX*NY*NZ;
        FluidField {
            energy: vec![epsilon_init; size],
            vx: vec![0.0; size],
            vy: vec![0.0; size],
            vz: vec![0.0; size],
            photon_density: vec![photon_init; size],
        }
    }
    
    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + NX*(y + NY*z)
    }

    fn laplacian(&self, arr: &Vec<f64>, x: usize, y: usize, z: usize) -> f64 {
        let val = arr[self.idx(x,y,z)];
        let xp = arr[self.idx((x+1)%NX, y, z)];
        let xm = arr[self.idx((x+NX-1)%NX, y, z)];
        let yp = arr[self.idx(x, (y+1)%NY, z)];
        let ym = arr[self.idx(x, (y+NY-1)%NY, z)];
        let zp = arr[self.idx(x, y, (z+1)%NZ)];
        let zm = arr[self.idx(x, y, (z+NZ-1)%NZ)];
        (xp+xm+yp+ym+zp+zm - 6.0*val)/(DX*DX)
    }
}

// A simple function simulating a "gravitational wave" metric factor variation
fn metric_factor(t: f64, x: f64) -> f64 {
    // A small oscillation in the metric
    // g_ii = a^2(t)*(1 + h cos(omega t - kx))
    let h = 1e-6; // tiny amplitude
    let omega = 2.0*PI; // 1 Hz wave
    let k = 2.0*PI/10.0; // wavelength 10 m
    1.0 + h*(omega*t - k*x).cos()
}

fn main() {
    let epsilon_init = 1e5; // J/m^3, arbitrary energy density
    let photon_init = 1e10; // photons/m^3, arbitrary
    let mut field = FluidField::new(epsilon_init, photon_init);

    let mut file = std::fs::File::create("fluid_lattice.csv").unwrap();
    writeln!(file, "time(s),avg_energy(J/m^3),avg_photon(m^-3)").unwrap();

    let steps = (TOTAL_TIME/DT) as usize;

    for step in 0..steps {
        let t = step as f64 * DT;
        // Evolve the fluid:
        let mut new_energy = field.energy.clone();
        let mut new_photons = field.photon_density.clone();

        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx = field.idx(x,y,z);
                    let e = field.energy[idx];
                    let p = qcd_pressure(e);

                    // Simple diffusion-like update for photons:
                    let lap_ph = field.laplacian(&field.photon_density, x,y,z);
                    let D_ph = 0.1; // chosen small diffusion coefficient
                    let n_ph = field.photon_density[idx];

                    // Include metric factor:
                    let alpha = metric_factor(t, x as f64 * DX);

                    // Update photon density (as if slightly affected by metric change)
                    // dn_ph/dt ~ D_ph * lap(n_ph)*alpha
                    let dn_ph = D_ph * lap_ph * DT * alpha;

                    // Energy density might be slightly modulated:
                    // In reality, you'd solve full fluid eq. Here we do a toy model:
                    let lap_e = field.laplacian(&field.energy, x,y,z);
                    let D_e = 0.01;
                    let de = D_e*lap_e*DT*alpha;

                    new_photons[idx] = n_ph + dn_ph;
                    new_energy[idx] = e + de - p*1e-3*DT; // a trivial source/sink term to mimic expansion
                }
            }
        }

        field.energy = new_energy;
        field.photon_density = new_photons;

        let avg_e = field.energy.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_ph = field.photon_density.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        writeln!(file, "{},{},{}", t, avg_e, avg_ph).unwrap();
    }

    println!("Fluid lattice simulation completed. Results in fluid_lattice.csv");
}
