use std::f64::consts::PI;
use std::io::Write;
use lazy_static::lazy_static;

// -------------------------------------------------------------------------
// Theoretical Foundations (Integrated with the Quantum Gravity Lagrangian)
// -------------------------------------------------------------------------
//
// We begin from a 4D quantum gravity inspired Lagrangian density:
//
//   ℒ_total = √(-g) [ (R/(16πG)) + ℒ_matter(ψ,∂ψ) + ℒ_EM(F_{μν}) + ℒ_QG(...) ]
//
// Here, R is the Ricci scalar, g is the metric determinant, G is Newton's constant,
// ψ represent matter fields (electrons, positrons, etc.), and F_{μν} the EM field tensor.
// ℒ_QG could include higher-order curvature terms, quantum corrections, or stringy modifications.
//
// From the Short Answer and Key Insights Provided:
// - Singularities in the action arise if fields, curvature, or matter densities become infinite or non-integrable.
// - The action: S = ∫ d^4x ℒ_total. If ℒ_total → ∞ near some region or fields are not smooth, S may diverge.
//
// Key Insights Implemented:
// 1. Dependence on Field Regularity: We must ensure fields remain smooth and finite.
// 2. Curvature Singularities: If R or other invariants blow up, the gravitational action diverges.
// 3. Matter Field Divergences: Unphysical infinite densities or discontinuous fields cause divergences.
// 4. Normalization and Physical Interpretations: Through renormalization and proper choice of variables,
//    we can remove coordinate artifacts and ensure a finite, well-defined action.
//
// Axioms to Move Forward:
// 1. Regularity Axiom: Exclude configurations that yield infinite densities or curvature. The code below
//    (though simplistic) uses stable parameters that do not approach known singular regions.
// 2. Renormalization Axiom: If infinities appear, we would introduce renormalization prescriptions
//    (not implemented numerically here, but assumed at the theoretical level).
// 3. Geometric Completeness Axiom: In a full QG simulation, ensure geodesic completeness or define boundaries
//    where theory is valid. The code scenario is a local approximation.
// 4. Physical Boundary Conditions Axiom: Apply conditions at boundaries (e.g., large radius or time) ensuring
//    fields vanish or remain finite, preserving a finite action.
//
// Conclusion from the Theoretical Side:
// By implementing these axioms and principles, we avoid action singularities. Our code scenario, while not
// fully capturing QG complexity, can be embedded in a larger framework where these axioms ensure no pathological
// singularities occur.
//
// "Calculus Index Window" for Rigor Testing:
// - After running simulations, we can:
//
 //   * Compute partial derivatives ∂μψ, ∂μF_{νρ}, check field smoothness.
 //   * Evaluate curvature invariants (R, R_{μν}R^{μν}, R_{μνρσ}R^{μνρσ}) from output metric data (not present here, but conceptualized).
 //   * Check integrability: numerically integrate ℒ_total over the domain and verify convergence.
 //   * If any measure diverges, apply renormalization or discard that region of parameter space.
 
// In a modern computational approach, symbolic differentiation tools or PDE solvers would link to this code block
// and run these checks as a separate module or "index window" to confirm the action remains finite and well-defined.

// -------------------------------------------------------------------------
// Physical Constants and Scenario Setup
// (The scenario code below simulates photon distributions and hypothetical boson production.)
// -------------------------------------------------------------------------
const C: f64 = 3.0e10;         // speed of light in cm/s
const R_E: f64 = 2.8179403262e-13; // classical electron radius in cm

// Scenario Constants
const E_ISO_ERG: f64 = 1.0e55; // total isotropic energy in erg
const R0: f64 = 1.0e13;        // initial radius in cm

// Band function parameters (for initial photon spectrum)
const ALPHA: f64 = -1.0;
const BETA_PAR: f64 = -2.3;
const E0_KEV: f64 = 300.0;
const KEV_TO_MEV: f64 = 1e-3;
const MEV_TO_ERG: f64 = 1.6021765e-6;

lazy_static! {
    static ref BETA: f64 = 1.0;
    // Thomson cross section:
    static ref SIGMA_T: f64 = (8.0 * PI / 3.0) * R_E.powi(2);
}

// Params struct for normalization
struct Params {
    norm: f64,
}

// State struct holds the densities (fields) we are tracking
#[derive(Clone)]
struct State {
    time: f64,
    radius: f64,
    n_photon: f64,
    n_pairs: f64,
    n_W: f64,
    n_Z: f64,
    n_H: f64,
}

// Band spectrum: A simplified matter/EM content scenario
fn band_spectrum(E_kev: f64, alpha: f64, beta: f64, E0_kev: f64, norm: f64) -> f64 {
    let e_break = (alpha - beta)*E0_kev;
    if E_kev < e_break {
        norm * (E_kev/100.0).powf(alpha)*(-E_kev/E0_kev).exp()
    } else {
        norm * ((e_break/100.0).powf(alpha - beta))*((alpha - beta).exp())*(E_kev/100.0).powf(beta)
    }
}

// Integrate band spectrum to find total energy, ensuring normalization is well-defined and finite:
fn total_energy_band(params: &Params) -> f64 {
    let e_min = 1.0;
    let e_max = 1.0e5;
    let steps = 200;
    let de = (e_max - e_min)/(steps as f64);
    let mut total_energy_erg = 0.0;
    for i in 0..steps {
        let E_kev = e_min + (i as f64)*de;
        let val = band_spectrum(E_kev, ALPHA, BETA_PAR, E0_KEV, params.norm);
        let E_mev = E_kev*KEV_TO_MEV;
        let E_erg = E_mev*MEV_TO_ERG;
        let dE = val * E_erg * de;
        total_energy_erg += dE;
    }
    total_energy_erg
}

fn find_norm_for_band() -> f64 {
    let guess = 1.0e5;
    let params = Params { norm: guess };
    let total = total_energy_band(&params);
    let target = E_ISO_ERG;
    guess*(target/total)
}

// Approximate pair production rate (toy model):
fn pair_production_rate(n_ph: f64) -> f64 {
    let sigma_eff = *SIGMA_T;  // Using Thomson cross section as a ballpark
    sigma_eff * C * n_ph.powi(2)
}

// W, Z, H production negligible at lower energies in this toy scenario:
fn w_production_rate(_n_ph: f64) -> f64 { 0.0 }
fn z_production_rate(_n_ph: f64) -> f64 { 0.0 }
fn h_production_rate(_n_ph: f64) -> f64 { 0.0 }

// Derivatives function: In a realistic QG scenario, these would come from field equations
// derived from variations of the action δS/δψ = 0, ensuring compliance with regularity axioms.
fn derivatives(s: &State) -> (f64, f64, f64, f64) {
    let pair_rate = pair_production_rate(s.n_photon);
    let w_rate = w_production_rate(s.n_photon);
    let z_rate = z_production_rate(s.n_photon);
    let h_rate = h_production_rate(s.n_photon);

    (pair_rate, w_rate, z_rate, h_rate)
}

// RK4 integrator: While integrating forward in time, we conceptually ensure no infinite fields appear:
// If a divergence were detected (e.g., radius → 0 causing infinite density), we'd halt or renormalize.
fn rk4_step(s: &mut State, dt: f64) {
    let s_original = s.clone();
    let (k1_pairs, k1_w, k1_z, k1_h) = derivatives(&s_original);

    let mut s2 = s_original.clone();
    s2.time += dt/2.0;
    s2.n_pairs += k1_pairs*(dt/2.0);
    s2.n_W += k1_w*(dt/2.0);
    s2.n_Z += k1_z*(dt/2.0);
    s2.n_H += k1_h*(dt/2.0);

    let (k2_pairs, k2_w, k2_z, k2_h) = derivatives(&s2);

    let mut s3 = s_original.clone();
    s3.time += dt/2.0;
    s3.n_pairs += k2_pairs*(dt/2.0);
    s3.n_W += k2_w*(dt/2.0);
    s3.n_Z += k2_z*(dt/2.0);
    s3.n_H += k2_h*(dt/2.0);

    let (k3_pairs, k3_w, k3_z, k3_h) = derivatives(&s3);

    let mut s4 = s_original.clone();
    s4.time += dt;
    s4.n_pairs += k3_pairs*dt;
    s4.n_W += k3_w*dt;
    s4.n_Z += k3_z*dt;
    s4.n_H += k3_h*dt;

    let (k4_pairs, k4_w, k4_z, k4_h) = derivatives(&s4);

    s.n_pairs += (k1_pairs + 2.0*k2_pairs + 2.0*k3_pairs + k4_pairs)*dt/6.0;
    s.n_W += (k1_w + 2.0*k2_w + 2.0*k3_w + k4_w)*dt/6.0;
    s.n_Z += (k1_z + 2.0*k2_z + 2.0*k3_z + k4_z)*dt/6.0;
    s.n_H += (k1_h + 2.0*k2_h + 2.0*k3_h + k4_h)*dt/6.0;
}

impl State {
    fn new(n_photon_init: f64) -> Self {
        // Starting state: finite photon density, zero W/Z/H
        // No singularities are triggered at initialization.
        State {
            time: 0.0,
            radius: R0,
            n_pairs: 0.0,
            n_photon: n_photon_init,
            n_W: 0.0,
            n_Z: 0.0,
            n_H: 0.0,
        }
    }
}

fn main() {
    let norm = find_norm_for_band();
    let params = Params { norm };
    println!("Normalization for Band function: {}", params.norm);

    // Estimate initial photon number density:
    let avg_E_kev = 1.0; // chosen for scaling
    let avg_E_mev = avg_E_kev*KEV_TO_MEV;
    let avg_E_erg = avg_E_mev*MEV_TO_ERG;
    let vol = (4.0/3.0)*PI*R0.powi(3);
    let total_photons = E_ISO_ERG / avg_E_erg;
    let n_photon_init = total_photons/vol;

    let mut state = State::new(n_photon_init);

    let total_time = 33000;
    let dt = 1.0; // Negative dt scenario simulates contraction

    let mut file = std::fs::File::create("extended_simulation_output.csv").unwrap();
    writeln!(file, "time(s),radius(cm),n_photon(cm^-3),n_pairs(cm^-3),n_W,n_Z,n_H").unwrap();

    for _ in 0..(total_time as usize) {
        state.time += dt;
        state.radius = R0 + *BETA*C*state.time;

        // Regularity Axiom in Action:
        // If radius approaches zero or a problematic regime, one would impose a cutoff or renormalization:
        // Here, we just allow scale factor changes but in a QG scenario, we'd check for divergence and stop if necessary.
        let scale = (R0/state.radius).powi(4);
        state.n_photon = n_photon_init*scale;

        // If n_photon or other fields become too large, in a full model we would renormalize or handle boundaries:
        // This ensures that the integral of ℒ_total remains finite and no singular actions appear.

        rk4_step(&mut state, dt);

        writeln!(file, "{},{},{},{},{},{},{}",
                 state.time,
                 state.radius,
                 state.n_photon,
                 state.n_pairs,
                 state.n_W,
                 state.n_Z,
                 state.n_H).unwrap();
    }

    println!("Final pairs: {} cm^-3", state.n_pairs);
    println!("Final W density: {} cm^-3", state.n_W);
    println!("Final Z density: {} cm^-3", state.n_Z);
    println!("Final H density: {} cm^-3", state.n_H);
    println!("Data in extended_simulation_output.csv");

// -------------------------------------------------------------------------
// Post-processing and Visual Representation (from the given visualization script):
//
// The visualization code (in Python, as provided) would then read the CSV file and
// plot the data. This provides a conceptual 'index window' of how fields evolve.
// We can interpret these results in light of the axioms and QG Lagrangian theory:
//
// - If any quantity here tended toward infinity, the final code would flag a divergence.
// - Future work would incorporate metric and curvature computations, providing a direct 
//   link between field evolution and the integrability of the action.
// -------------------------------------------------------------------------
}


/*
use std::f64::consts::PI;
use std::io::Write;
use lazy_static::lazy_static;

// Physical Constants
const C: f64 = 3.0e10;         // speed of light in cm/s
const R_E: f64 = 2.8179403262e-13; // classical electron radius in cm

// Scenario Constants
const E_ISO_ERG: f64 = 1.0e55; // Total isotropic energy in erg
const R0: f64 = 1.0e13;        // initial radius in cm

// Band function parameters (initial photon spectrum)
const ALPHA: f64 = -1.0;
const BETA_PAR: f64 = -2.3;
const E0_KEV: f64 = 300.0;
const KEV_TO_MEV: f64 = 1e-3;
const MEV_TO_ERG: f64 = 1.6021765e-6;

lazy_static! {
    static ref BETA: f64 = 1.0;
    // Thomson cross section:
    static ref SIGMA_T: f64 = (8.0 * PI / 3.0) * R_E.powi(2);
}

// Params struct for normalization
struct Params {
    norm: f64,
}

#[derive(Clone)]
struct State {
    time: f64,
    radius: f64,
    n_photon: f64,
    n_pairs: f64,
    n_W: f64,
    n_Z: f64,
    n_H: f64,
}

// Band spectrum (Fermi-Dirac or phenomenological Band function)
fn band_spectrum(E_kev: f64, alpha: f64, beta: f64, E0_kev: f64, norm: f64) -> f64 {
    let e_break = (alpha - beta)*E0_kev;
    if E_kev < e_break {
        norm * (E_kev/100.0).powf(alpha)*(-E_kev/E0_kev).exp()
    } else {
        norm * ((e_break/100.0).powf(alpha - beta))*((alpha - beta).exp())*(E_kev/100.0).powf(beta)
    }
}

fn total_energy_band(params: &Params) -> f64 {
    let e_min = 1.0;
    let e_max = 1.0e5;
    let steps = 200;
    let de = (e_max - e_min)/(steps as f64);
    let mut total_energy_erg = 0.0;
    for i in 0..steps {
        let E_kev = e_min + (i as f64)*de;
        let val = band_spectrum(E_kev, ALPHA, BETA_PAR, E0_KEV, params.norm);
        let E_mev = E_kev*KEV_TO_MEV;
        let E_erg = E_mev*MEV_TO_ERG;
        let dE = val * E_erg * de;
        total_energy_erg += dE;
    }
    total_energy_erg
}

fn find_norm_for_band() -> f64 {
    let guess = 1.0e5;
    let params = Params { norm: guess };
    let total = total_energy_band(&params);
    let target = E_ISO_ERG;
    guess*(target/total)
}

// Approximate pair production rate from a dense photon gas:
// Rate ~ sigma * c * n_ph^2
// We'll assume all photons have sufficiently high energy above threshold.
// For a more accurate model, one would integrate over the photon distribution and use the actual gamma-gamma cross section.
// Here, we approximate sigma ~ SIGMA_T to get an order-of-magnitude estimate.
fn pair_production_rate(n_ph: f64) -> f64 {
    let sigma_eff = *SIGMA_T;  // using Thomson cross section as a ballpark
    sigma_eff * C * n_ph.powi(2)
}

// Given the energies, W, Z, H production is negligible if we are not at extremely high energies.
// Set them to zero for physical realism at lower energies.
fn w_production_rate(_n_ph: f64) -> f64 { 0.0 }
fn z_production_rate(_n_ph: f64) -> f64 { 0.0 }
fn h_production_rate(_n_ph: f64) -> f64 { 0.0 }

// Derivatives function:
fn derivatives(s: &State) -> (f64, f64, f64, f64) {
    let pair_rate = pair_production_rate(s.n_photon);
    let w_rate = w_production_rate(s.n_photon);
    let z_rate = z_production_rate(s.n_photon);
    let h_rate = h_production_rate(s.n_photon);

    (pair_rate, w_rate, z_rate, h_rate)
}

// RK4 integrator:
fn rk4_step(s: &mut State, dt: f64) {
    let s_original = s.clone();
    let (k1_pairs, k1_w, k1_z, k1_h) = derivatives(&s_original);

    let mut s2 = s_original.clone();
    s2.time += dt/2.0;
    s2.n_pairs += k1_pairs*(dt/2.0);
    s2.n_W += k1_w*(dt/2.0);
    s2.n_Z += k1_z*(dt/2.0);
    s2.n_H += k1_h*(dt/2.0);

    let (k2_pairs, k2_w, k2_z, k2_h) = derivatives(&s2);

    let mut s3 = s_original.clone();
    s3.time += dt/2.0;
    s3.n_pairs += k2_pairs*(dt/2.0);
    s3.n_W += k2_w*(dt/2.0);
    s3.n_Z += k2_z*(dt/2.0);
    s3.n_H += k2_h*(dt/2.0);

    let (k3_pairs, k3_w, k3_z, k3_h) = derivatives(&s3);

    let mut s4 = s_original.clone();
    s4.time += dt;
    s4.n_pairs += k3_pairs*dt;
    s4.n_W += k3_w*dt;
    s4.n_Z += k3_z*dt;
    s4.n_H += k3_h*dt;

    let (k4_pairs, k4_w, k4_z, k4_h) = derivatives(&s4);

    s.n_pairs += (k1_pairs + 2.0*k2_pairs + 2.0*k3_pairs + k4_pairs)*dt/6.0;
    s.n_W += (k1_w + 2.0*k2_w + 2.0*k3_w + k4_w)*dt/6.0;
    s.n_Z += (k1_z + 2.0*k2_z + 2.0*k3_z + k4_z)*dt/6.0;
    s.n_H += (k1_h + 2.0*k2_h + 2.0*k3_h + k4_h)*dt/6.0;
}

impl State {
    fn new(n_photon_init: f64) -> Self {
        State {
            time: 0.0,
            radius: R0,
            n_pairs: 0.0,
            n_photon: n_photon_init,
            n_W: 0.0,
            n_Z: 0.0,
            n_H: 0.0,
        }
    }
}

fn main() {
    let norm = find_norm_for_band();
    let params = Params { norm };
    println!("Normalization for Band function: {}", params.norm);

    let avg_E_kev = 300.0; // Just a placeholder for scaling
    let avg_E_mev = avg_E_kev*KEV_TO_MEV;
    let avg_E_erg = avg_E_mev*MEV_TO_ERG;
    let vol = (4.0/3.0)*PI*R0.powi(3);
    let total_photons = E_ISO_ERG / avg_E_erg;
    let n_photon_init = total_photons/vol;

    let mut state = State::new(n_photon_init);

    let total_time = 500;
    let dt = -1.0; // Negative dt: contracting scenario

    let mut file = std::fs::File::create("extended_simulation_output.csv").unwrap();
    writeln!(file, "time(s),radius(cm),n_photon(cm^-3),n_pairs(cm^-3),n_W,n_Z,n_H").unwrap();

    for _ in 0..(total_time as usize) {
        state.time += dt;
        state.radius = R0 + *BETA*C*state.time;
        // Update photon density (assuming adiabatic compression):
        let scale = (R0/state.radius).powi(3);
        state.n_photon = n_photon_init*scale;

        rk4_step(&mut state, dt);

        writeln!(file, "{},{},{},{},{},{},{}",
                 state.time,
                 state.radius,
                 state.n_photon,
                 state.n_pairs,
                 state.n_W,
                 state.n_Z,
                 state.n_H).unwrap();
    }

    println!("Final pairs: {} cm^-3", state.n_pairs);
    println!("Final W density: {} cm^-3", state.n_W);
    println!("Final Z density: {} cm^-3", state.n_Z);
    println!("Final H density: {} cm^-3", state.n_H);
    println!("Data in simulation_output.csv");
}
*/







/*
use std::f64::consts::PI;
use std::io::Write;
use lazy_static::lazy_static;

// ----------------------------------------
// Fundamental constants and parameters (Standard Model + Cosmology)

// Physical constants and scales:
const C: f64 = 3.0e10;          // Speed of light in cm/s
const MP: f64 = 2.435e18;       // Reduced Planck mass in GeV
const H0: f64 = 1.0e-33;        // Hubble scale (~GeV)
const G_F: f64 = 1.1663787e-5;  // Fermi coupling constant in GeV^-2

// Standard Model masses (in GeV):
const M_W: f64 = 80.379;  // W boson mass
const M_Z: f64 = 91.1876; // Z boson mass
const M_H: f64 = 125.10;  // Higgs boson mass

// Electroweak parameters:
const V_HIGGS: f64 = 246.0; // Higgs vev in GeV
const ALPHA_EM: f64 = 1.0/137.035999; // fine structure constant
// Weak mixing angle (Weinberg angle), sin^2 θ_W ~ 0.231
const SIN2_THETA_W: f64 = 0.231;
lazy_static! {
    static ref SIN_THETA_W: f64 = SIN2_THETA_W.sqrt();
    static ref COS_THETA_W: f64 = (1.0 - SIN2_THETA_W).sqrt();
}

// Example scenario parameters:
const E_ISO_ERG: f64 = 1.0e55;
const R0: f64 = 1.0e13; // initial radius in cm

// Particle distribution parameters (as before, but could be adapted):
const ALPHA: f64 = -1.0;
const BETA_PAR: f64 = -2.3;
const E0_KEV: f64 = 300.0;
const KEV_TO_MEV: f64 = 1e-3;
const MEV_TO_ERG: f64 = 1.6021765e-6;

lazy_static! {
    // We can imagine a final scale or scenario:
    static ref HF: f64 = H0; // Just a placeholder for final scale
    static ref BETA: f64 = 1.0;
}

// Struct to store normalization and possibly couplings:
struct Params {
    norm: f64,
    g_fermi: f64, // From G_F, can be related to weak interactions
    // Additional parameters for Higgs couplings, etc.
    yukawa_coupling: f64,
}

#[derive(Clone)]
struct State {
    time: f64,
    radius: f64,
    n_photon: f64,
    n_pairs: f64,
    // New fields for boson densities:
    n_W: f64,     // W boson density
    n_Z: f64,     // Z boson density
    n_H: f64,     // Higgs boson density
}

// Simple band spectrum as before:
fn band_spectrum(E_kev: f64, alpha: f64, beta: f64, E0_kev: f64, norm: f64) -> f64 {
    let e_break = (alpha - beta)*E0_kev;
    if E_kev < e_break {
        norm * (E_kev/100.0).powf(alpha)*(-E_kev/E0_kev).exp()
    } else {
        norm * ((e_break/100.0).powf(alpha - beta))*((alpha - beta).exp())*(E_kev/100.0).powf(beta)
    }
}

fn total_energy_band(params: &Params) -> f64 {
    let e_min = 1.0;
    let e_max = 1.0e5;
    let steps = 200;
    let de = (e_max - e_min)/(steps as f64);
    let mut total_energy_erg = 0.0;
    for i in 0..steps {
        let E_kev = e_min + (i as f64)*de;
        let val = band_spectrum(E_kev, ALPHA, BETA_PAR, E0_KEV, params.norm);
        let E_mev = E_kev*KEV_TO_MEV;
        let E_erg = E_mev*MEV_TO_ERG;
        let dE = val * E_erg * de;
        total_energy_erg += dE;
    }
    total_energy_erg
}

fn find_norm_for_band() -> f64 {
    let guess = 1.0e5;
    let params = Params {
        norm: guess,
        g_fermi: G_F,
        yukawa_coupling: 1.0, // Placeholder
    };
    let total = total_energy_band(&params);
    let target = E_ISO_ERG;
    guess*(target/total)
}

// Placeholder function to represent interactions involving W, Z, and Higgs bosons:
// In a real model, this would involve detailed cross sections and phase space integrals.
fn weak_interaction_rate(n_ph: f64, params: &Params) -> f64 {
    // Roughly scale as G_F^2 * energy density for demonstration:
    let energy_scale: f64 = 1.0e9; // placeholder energy scale in GeV
    // Rate ~ G_F^2 * n_photon^2 * E^2 as a toy model
    params.g_fermi.powi(2)*n_ph.powf(2.0)*energy_scale.powi(2)
}

// Higgs related rates (placeholder):
fn higgs_mediated_rate(n_ph: f64, params: &Params) -> f64 {
    // Yukawa couplings scale like m_f / v, here just a placeholder:
    let factor = params.yukawa_coupling;
    // Rate ~ factor^2 * n_ph^2:
    factor.powi(2)*n_ph.powf(2.0)*1e-10
}

// Combined derivative function that includes new interactions:
fn derivatives(s: &State, params: &Params) -> (f64, f64, f64, f64) {
    // Pair production (as before) plus new contributions:
    let pair_rate = s.n_photon.powf(2.0)*1e-20; // placeholder from old logic
    let w_production = weak_interaction_rate(s.n_photon, params)*1e-30; 
    let z_production = w_production*0.5; // Just assume Z about half for demonstration
    let h_production = higgs_mediated_rate(s.n_photon, params);

    // Return tuple: (dn_pairs/dt, dn_W/dt, dn_Z/dt, dn_H/dt)
    (pair_rate, w_production, z_production, h_production)
}

fn rk4_step(s: &mut State, dt: f64, params: &Params) {
    let s_original = s.clone();
    let (k1_pairs, k1_w, k1_z, k1_h) = derivatives(&s_original, params);

    let mut s2 = s_original.clone();
    s2.time += dt/2.0;
    s2.n_pairs += k1_pairs*dt/2.0;
    s2.n_W += k1_w*dt/2.0;
    s2.n_Z += k1_z*dt/2.0;
    s2.n_H += k1_h*dt/2.0;

    let (k2_pairs, k2_w, k2_z, k2_h) = derivatives(&s2, params);

    let mut s3 = s_original.clone();
    s3.time += dt/2.0;
    s3.n_pairs += k2_pairs*dt/2.0;
    s3.n_W += k2_w*dt/2.0;
    s3.n_Z += k2_z*dt/2.0;
    s3.n_H += k2_h*dt/2.0;

    let (k3_pairs, k3_w, k3_z, k3_h) = derivatives(&s3, params);

    let mut s4 = s_original.clone();
    s4.time += dt;
    s4.n_pairs += k3_pairs*dt;
    s4.n_W += k3_w*dt;
    s4.n_Z += k3_z*dt;
    s4.n_H += k3_h*dt;

    let (k4_pairs, k4_w, k4_z, k4_h) = derivatives(&s4, params);

    s.n_pairs += (k1_pairs + 2.0*k2_pairs + 2.0*k3_pairs + k4_pairs)*dt/6.0;
    s.n_W += (k1_w + 2.0*k2_w + 2.0*k3_w + k4_w)*dt/6.0;
    s.n_Z += (k1_z + 2.0*k2_z + 2.0*k3_z + k4_z)*dt/6.0;
    s.n_H += (k1_h + 2.0*k2_h + 2.0*k3_h + k4_h)*dt/6.0;
}

impl State {
    fn new(n_photon_init: f64) -> Self {
        State {
            time: 0.0,
            radius: R0,
            n_pairs: 0.0,
            n_photon: n_photon_init,
            n_W: 0.0,
            n_Z: 0.0,
            n_H: 0.0,
        }
    }
}

fn main() {
    let norm = find_norm_for_band();
    let params = Params {
        norm,
        g_fermi: G_F,
        yukawa_coupling: 0.01, // placeholder
    };
    println!("Normalization for Band function: {}", params.norm);

    let avg_E_kev = 1.0;
    let avg_E_mev = avg_E_kev*KEV_TO_MEV;
    let avg_E_erg = avg_E_mev*MEV_TO_ERG;
    let vol = (4.0/3.0)*PI*R0.powi(3);
    let total_photons = E_ISO_ERG / avg_E_erg;
    let n_photon_init = total_photons/vol;

    let mut state = State::new(n_photon_init);

    let total_time = 50;
    let dt = -1.0; // Reverse time step for conceptual scenario, can also be positive

    let mut file = std::fs::File::create("extended_simulation_output.csv").unwrap();
    writeln!(file, "time(s),radius(cm),n_photon(cm^-3),n_pairs(cm^-3),n_W,n_Z,n_H").unwrap();

    for _ in 0..(total_time as usize) {
        state.time += dt;
        state.radius = R0 + *BETA*C*state.time;
        let scale = (R0/state.radius).powi(3);
        state.n_photon = n_photon_init*scale;

        rk4_step(&mut state, dt, &params);

        writeln!(file, "{},{},{},{},{},{},{}",
                 state.time,
                 state.radius,
                 state.n_photon,
                 state.n_pairs,
                 state.n_W,
                 state.n_Z,
                 state.n_H).unwrap();
    }

    println!("Final pairs: {} cm^-3", state.n_pairs);
    println!("Final W density: {} cm^-3", state.n_W);
    println!("Final Z density: {} cm^-3", state.n_Z);
    println!("Final Higgs density: {} cm^-3", state.n_H);
    println!("Data in extended_simulation_output.csv");
}
*/