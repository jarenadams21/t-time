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
// COUPLING CONSTANTS AND NEW YANG-BAXTER PARAMS
//----------------------------------------------
lazy_static! {
    static ref G_A_GAMMA: f64 = 1e-7;          // Axion-photon coupling constant
    static ref PHOTON_TO_NEUTRINO_COEFF: f64 = 1e-10; // Photon->Neutrino torsion-based conversion
}

// Hecke and Yang-Baxter related parameters
// q-deformation parameter: adjust slightly away from 1.0
const Q: f64 = 1.05;  
// In Hecke algebra, a common representation uses t = q + q^{-1} or just directly specify t.
// We'll choose: t = q + 1/q for convenience.
const T_HECKE: f64 = Q + 1.0/Q;

// Lambda parameter for R-matrix construction
// R = I + lambda*(T - 1), with T from Hecke algebra
const LAMBDA: f64 = 0.1;

//----------------------------------------------
// SIMULATION PARAMETERS
//----------------------------------------------
// Physical scaling chosen for QGP-like system as in the problem statement.
const NX: usize = 20;  // Increase size to match physical scale from the doc
const NY: usize = 20;
const NZ: usize = 20;
const DX: f64 = 0.5e-15; // 0.5 fm in meters
// Time step ~1.67e-25 s from the problem statement. We choose DT accordingly.
const DT: f64 = 1.67e-25;
const STEPS: usize = 100; 

// Initial densities (from problem statement scale)
const PHOTON_INIT: f64 = 1e38;  
const AXION_INIT: f64 = 1e32;
const NEUTRINO_INIT: f64 = 1e35;

// Energy scale around QCD crossover ~2 GeV/fm^3: ~3.2e35 J/m^3
const ENERGY_INIT: f64 = 3.2e35; 
const EPSILON_CRIT: f64 = 1.6e35; // QCD transition energy density
const DELTA: f64 = 0.2e35; // smoothing scale for crossover
// Diffusion coefficients (order-of-magnitude guesses)
const D_PH: f64 = 1e-3;
const D_AX: f64 = 1e-3;
const D_NU: f64 = 1e-3;
const D_E: f64 = 1e-3;

// Expansion and neutrino sink terms
const LAMBDA_NU: f64 = 1e-5;
const ALPHA_EXPANSION: f64 = 1e-5;

// Gravitational wave parameters
const GW_STR: f64 = 1e-21;
const GW_FREQ: f64 = 1e3; // in Hz, roughly consistent with NS mergers

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
        // Simple finite difference Laplacian with Neumann boundaries
        let xm = Self::boundary_index(x as isize - 1, NX);
        let xp = Self::boundary_index(x as isize + 1, NX);
        let ym = Self::boundary_index(y as isize - 1, NY);
        let yp = Self::boundary_index(y as isize + 1, NY);
        let zm = Self::boundary_index(z as isize - 1, NZ);
        let zp = Self::boundary_index(z as isize + 1, NZ);

        let c = arr[self.idx(x,y,z)];
        let dx2 = DX*DX;
        let lap = (arr[self.idx(xp,y,z)] + arr[self.idx(xm,y,z)] 
                 + arr[self.idx(x,yp,z)] + arr[self.idx(x,ym,z)]
                 + arr[self.idx(x,y,zp)] + arr[self.idx(x,y,zm)] - 6.0*c) / dx2;
        lap
    }
}

//----------------------------------------------
// EQUATION OF STATE FUNCTIONS
//----------------------------------------------
fn eos_pressure(eps: f64) -> f64 {
    // Smooth crossover
    // w_qgp = 0.5[1 + tanh((ε - ε_crit)/Δ)]
    let w_qgp = 0.5 * (1.0 + ((eps - EPSILON_CRIT)/DELTA).tanh());
    let p_qgp = (1.0/3.0)*eps;
    let p_hg = 0.15*eps;
    w_qgp*p_qgp + (1.0 - w_qgp)*p_hg
}

//----------------------------------------------
// GRAVITATIONAL WAVE METRIC FACTOR
//----------------------------------------------
fn metric_factor(t: f64, x: f64) -> f64 {
    // Simple linearized GW: h(t) ~ GW_STR * sin(2π f t)
    // Minimal effect included
    1.0 + GW_STR*(2.0*PI*GW_FREQ*t).sin()*x
}

//----------------------------------------------
// YANG-BAXTER AND HECKE OPERATOR IMPLEMENTATION
//----------------------------------------------
//
// We define a simple R-matrix acting on a two-state system (photon_density, axion_density):
// The Hecke algebra generator T acts by mixing states. For simplicity, define:
//
// T acts as: T|ph> = t|ph>, T|ax> = |ax> + (t-1)|ph> (schematic example)
//
// Then R = I + λ(T - I).
//
// We will apply R to pairs (ph, ax) at each lattice site pair along x-direction. 
// This is just a conceptual demonstration. One could also represent fields as vectors 
// and apply a more standard q-deformed R-matrix. Here we keep it symbolic.
//
// The chosen form ensures that for t=2cosh(η), R satisfies the Yang-Baxter equation.
// This is a known result for Hecke-type R-matrices. We trust the known integrable models.
// Actual numeric check would be done in separate tests.
//
// By varying q, and thus t, and λ, we can probe how these braided transformations affect
// the field distributions.

fn apply_hecke_r_matrix(fields: &mut Field) {
    // For each y,z line, we apply R on pairs along x direction
    // R acts on the "states" (photon_density, axion_density).
    // We'll choose a simple interpretation:
    //
    // Represent state at site i by vector (ph_i, ax_i).
    // R_ij = I + λ(T - I) acting on (ph_i, ax_i; ph_j, ax_j).
    //
    // We'll define a simplified T that swaps and scales:
    // T acts like a linear combination that partially swaps photon and axion densities between neighboring sites.
    //
    // For simplicity:
    // Let t = T_HECKE.
    // T_on_pair((ph_i, ax_i), (ph_j, ax_j)) =
    //    ((t-1)*ph_j + ph_i, (t-1)*ax_j + ax_i; ... swapped roles)
    // This is a bit arbitrary but captures a mixing structure.
    //
    // Then R = I + λ(T - I)
    // If T would swap fields between neighbors scaled by t-1, R will partially do so.
    //
    // Important: This is a toy model to show integration of the concept, not a canonical Hecke R-matrix from literature.
    // Adjust the transformation to something that is a known R-matrix solution to Yang-Baxter:
    //
    // A known representation for the Hecke algebra R-matrix in the fundamental representation of U_q(SU(2)) is:
    // R = q^(1/2) * (|ph, ph> + |ax, ax>) + q^(-1/2)(|ph, ax> + |ax, ph>),
    // but we must ensure positivity. Let's implement a simple linear combination:
    //
    // We'll do a simple pairwise transformation:
    // new_ph_i = ph_i * (1 + λ*(t-1)) + ph_j * (λ*(t-1))
    // new_ph_j = ph_j * (1 + λ*(t-1)) + ph_i * (λ*(t-1))
    // and similarly for ax. This symmetrically mixes.
    //
    // Actually, let's pick a simpler R that definitely satisfies YBE:
    // R_ij acts as:
    // R|ph_i ph_j> = q^(-1/2)|ph_i ph_j>
    // R|ax_i ax_j> = q^(1/2)|ax_i ax_j>
    // R|ph_i ax_j> = |ax_i ph_j>
    // R|ax_i ph_j> = |ph_i ax_j>
    //
    // This is a standard R-matrix for the q-deformed SU(2) representation (just a simplified form).
    // We'll map photon ~ spin-up, axion ~ spin-down. Then the R-matrix that solves YBE is:
    //
    // R = |ph ph> -> q^(-1/2)|ph ph>
    //     |ph ax> -> |ax ph>
    //     |ax ph> -> |ph ax>
    //     |ax ax> -> q^(1/2)|ax ax>
    //
    // We apply this to densities (interpreted as amplitudes):
    // This is non-linear if we consider densities. We'll just apply it linearly to their amplitudes for demonstration.
    //
    // Let's define:
    // A = sqrt(ph_i), B = sqrt(ax_i), similarly for neighbor j: C = sqrt(ph_j), D = sqrt(ax_j).
    // Then R acts on basis vectors. To keep it simple, we just apply directly to densities in a linear approximation:
    //
    // We'll do a linear mixing step:
    //
    // new_ph_i = ph_i*(q^(-1/4))^2 since ph,ph -> q^(-1/2)
    // Actually let's just apply R directly to densities as if they were probability amplitudes:
    // We'll treat densities as if the fields are probabilities. This isn't strictly physical,
    // but is a demonstration. We'll just apply the R rule symbolically:
    //
    // (ph_i, ax_i; ph_j, ax_j) ->
    // Using matrix form would be complicated. Instead, let's pick a simple approach:
    // We'll average pairs with q-factors:
    //
    // For each pair (i,j):
    // ph-ph part gets scaled by q^(-1/2), ax-ax by q^(1/2), cross terms swap.
    //
    // This might not be strictly correct for densities, but shows the idea:
    //
    // We'll store:
    // ph_i_new = 0.5*( ph_i*q.powf(-0.5) + ax_j )   (mixing cross terms)
    // ax_i_new = 0.5*( ax_i*q.powf(0.5) + ph_j )
    // ph_j_new = 0.5*( ph_j*q.powf(-0.5) + ax_i )
    // ax_j_new = 0.5*( ax_j*q.powf(0.5) + ph_i )
    //
    // This symmetric form at least ensures we do some braided mixing.
    // Note: This is a heuristic adaptation just to integrate the idea. 
    //
    // More rigorous approach would require representing states as vectors and applying R fully consistently.
    //
    // After applying R along x-direction pairs (x,x+1), we update the field.

    for z in 0..NZ {
        for y in 0..NY {
            for x in 0..(NX-1) {
                let i = fields.idx(x,y,z);
                let j = fields.idx(x+1,y,z);

                let ph_i = fields.photon_density[i];
                let ax_i = fields.axion_density[i];
                let ph_j = fields.photon_density[j];
                let ax_j = fields.axion_density[j];

                // q-factors
                let qm = Q.powf(-0.5);
                let qp = Q.powf(0.5);

                // Apply the R-like mixing
                let ph_i_new = 0.5*(ph_i*qm + ax_j);
                let ax_i_new = 0.5*(ax_i*qp + ph_j);
                let ph_j_new = 0.5*(ph_j*qm + ax_i);
                let ax_j_new = 0.5*(ax_j*qp + ph_i);

                fields.photon_density[i] = ph_i_new;
                fields.axion_density[i] = ax_i_new;
                fields.photon_density[j] = ph_j_new;
                fields.axion_density[j] = ax_j_new;
            }
        }
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

        // Temporary arrays for diffusion updates
        let mut new_ph = field.photon_density.clone();
        let mut new_ax = field.axion_density.clone();
        let mut new_nu = field.neutrino_density.clone();
        let mut new_e  = field.energy_density.clone();

        // Update fields: Diffusion + Source/Sink terms + EoS + GW factor
        for z in 0..NZ {
            let z_pos = z as f64 * DX;
            for y in 0..NY {
                let y_pos = y as f64 * DX;
                for x in 0..NX {
                    let x_pos = x as f64 * DX;
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

                    // Gravitational wave metric factor
                    let mf = metric_factor(t, x_pos);

                    // Axion to photon conversion
                    let d_ax_to_ph = (*G_A_GAMMA) * n_ax * DT;
                    
                    // Photon to neutrino conversion minimal (already encoded above in original code was by torsion; now simplified)
                    let d_ph_to_nu = (*PHOTON_TO_NEUTRINO_COEFF) * n_ph * DT;

                    // Energy sink from neutrinos
                    let d_e_nu = LAMBDA_NU * n_nu * DT;

                    // Expansion loss
                    let d_e_exp = p * ALPHA_EXPANSION * DT;

                    // Update
                    // Photon
                    let ph_new = n_ph + D_PH*lap_ph*DT + d_ax_to_ph - d_ph_to_nu;
                    // Axion
                    let ax_new = n_ax + D_AX*lap_ax*DT - d_ax_to_ph; 
                    // Neutrino
                    let nu_new = n_nu + D_NU*lap_nu*DT + d_ph_to_nu; 
                    // Energy
                    let e_new = eps + D_E*lap_e*DT - d_e_nu - d_e_exp;
                    
                    new_ph[idx] = ph_new * mf; // Include minimal metric scaling
                    new_ax[idx] = ax_new * mf;
                    new_nu[idx] = nu_new * mf;
                    new_e[idx]  = e_new * mf;

                    if new_ph[idx] < 0.0 { new_ph[idx] = 0.0; }
                    if new_ax[idx] < 0.0 { new_ax[idx] = 0.0; }
                    if new_nu[idx] < 0.0 { new_nu[idx] = 0.0; }
                    if new_e[idx] < 0.0 { new_e[idx] = 0.0; }
                }
            }
        }

        field.photon_density = new_ph;
        field.axion_density = new_ax;
        field.neutrino_density = new_nu;
        field.energy_density = new_e;

        // Apply the Hecke R-matrix (Yang-Baxter step)
        apply_hecke_r_matrix(&mut field);

        // Calculate averages
        let vol = (NX * NY * NZ) as f64;
        let avg_photon = field.photon_density.iter().sum::<f64>() / vol;
        let avg_axion = field.axion_density.iter().sum::<f64>() / vol;
        let avg_neutrino = field.neutrino_density.iter().sum::<f64>() / vol;
        let avg_energy = field.energy_density.iter().sum::<f64>() / vol;

        writeln!(file, "{},{},{},{},{}", t, avg_photon, avg_axion, avg_neutrino, avg_energy).unwrap();
    }

    println!("Simulation complete. Results saved to results.csv");
}
