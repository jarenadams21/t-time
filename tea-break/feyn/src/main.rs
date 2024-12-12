use std::f64::consts::PI;
use std::io::Write;
use lazy_static::lazy_static;

// ---------------------------------------------------------
// Physical Constants (CODATA and recognized standards ~2024)
// ---------------------------------------------------------
const C: f64 = 2.99792458e8;    // m/s
const H: f64 = 6.62607015e-34;  // J·s
const K_B: f64 = 1.380649e-23;  // J/K
const G: f64 = 6.67430e-11;     // m^3 kg^-1 s^-2

// Electron and charge-related constants (unused directly, but kept for completeness)
const M_E: f64 = 9.10938356e-31; // kg
const E_CHARGE: f64 = 1.602176634e-19; // C

// Classical electron radius (unused, but kept)
const R_E: f64 = 2.8179403262e-15; // m
lazy_static! {
    static ref SIGMA_T: f64 = (8.0*PI/3.0)*R_E.powi(2)*5.0/3.0; // Thomson cross-section modification (not directly used)
}

// We now focus on CMB parameters:
const T_CMB: f64 = 2.7255; // K ~ Current CMB temperature

// Earth's gravitational acceleration:
const G_EARTH: f64 = 9.81; // m/s^2

// Radius of spherical region of interest:
const R0: f64 = 10.0; // meters

// Compute photon number density for CMB:
// We integrate the Planck distribution for photon number density:
// n(ν)dν = (8 π ν² / c³) [1/(e^(hν/(k_B T)) - 1)] dν
fn planck_number_density(T: f64) -> f64 {
    let h = H;
    let kb = K_B;
    let c = C;
    let t = T;

    // Frequency range chosen to cover microwave domain up to infrared:
    // The CMB peaks at ~160 GHz (1.6e11 Hz), but we integrate a broad range:
    let nu_min = 1.0e7;   // 10 MHz, well below CMB peak
    let nu_max = 1.0e14;  // 100 THz, well above CMB tail
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

// Pair production negligible at these conditions:
fn pair_production_rate(_n_ph: f64) -> f64 {
    0.0
}

#[derive(Clone)]
struct State {
    time: f64,
    radius: f64,
    n_photon: f64,
    n_pairs: f64,
}

fn derivatives(s: &State) -> f64 {
    pair_production_rate(s.n_photon)
}

fn rk4_step(s: &mut State, dt: f64) {
    let s0 = s.clone();
    let k1 = derivatives(&s0);

    let mut s2 = s0.clone();
    s2.time += dt/2.0;
    s2.n_pairs += k1*(dt/2.0);
    let k2 = derivatives(&s2);

    let mut s3 = s0.clone();
    s3.time += dt/2.0;
    s3.n_pairs += k2*(dt/2.0);
    let k3 = derivatives(&s3);

    let mut s4 = s0.clone();
    s4.time += dt;
    s4.n_pairs += k3*dt;
    let k4 = derivatives(&s4);

    s.n_pairs += (k1 + 2.0*k2 + 2.0*k3 + k4)*dt/6.0;
}

impl State {
    fn new(n_photon_init: f64) -> Self {
        State {
            time: 0.0,
            radius: R0,
            n_photon: n_photon_init,
            n_pairs: 0.0,
        }
    }
}

// Photon-exclusive tunneling condition check function
// This is a placeholder function that may represent a condition derived from your theory.
// It's been left as-is, just changing the photon number density input to CMB-based results.
fn tunneling_condition(theta: f64, _phi: f64) -> f64 {
    let grav_correction = 1e-9 * G_EARTH;
    let ln_term = ((R0/10.0).powi(2) + 1.0).ln();
    let sqrt_ln = ln_term.sqrt();

    // Use the CMB photon number density:
    let n_ph = planck_number_density(T_CMB);
    let photon_correction = n_ph * 1e-36; // Scaled down to avoid overshadowing other terms

    let f = theta.cos() - sqrt_ln + photon_correction - grav_correction;
    f
}

fn main() {
    // Compute the photon number density for the CMB:
    let n_photon_init = planck_number_density(T_CMB);
    println!("CMB photon number density: {} photons/m^3", n_photon_init);

    // Initialize state and run a simulation for pair density (though we expect no pairs):
    let mut state = State::new(n_photon_init);

    // Integrate over about 730 seconds (~12 minutes) as a placeholder:
    let total_time_s = 7300000.0; 
    let dt = 1.0;
    let mut file = std::fs::File::create("cmb_simulation_output.csv").unwrap();
    writeln!(file, "time(s),radius(m),n_photon(m^-3),n_pairs(m^-3)").unwrap();

    for _ in 0..(total_time_s as usize) {
        rk4_step(&mut state, dt);
        state.time += dt;
        writeln!(file, "{},{},{},{}",
                 state.time,
                 state.radius,
                 state.n_photon,
                 state.n_pairs).unwrap();
    }

    println!("Final pair density: {} m^-3", state.n_pairs);
    println!("Data saved to cmb_simulation_output.csv");

    // Search over directions for tunneling condition under CMB conditions:
    let n_theta = 320; //98
    let n_phi = 450; //177
    let tolerance = 1e-6;
    let mut tunneling_directions = Vec::new();

    for i in 0..n_theta {
        for j in 0..n_phi {
            let theta = (i as f64)/(n_theta as f64)*PI;  // 0 to π
            let phi = (j as f64)/(n_phi as f64)*2.0*PI;  // 0 to 2π
            let val = tunneling_condition(theta, phi);
            if val.abs() < tolerance {
                tunneling_directions.push((theta, phi, val));
            }
        }
    }

    if tunneling_directions.is_empty() {
        println!("No photon-exclusive tunneling directions found under CMB conditions.");
    } else {
        println!("Photon-exclusive tunneling possible in {} directions under CMB conditions!", tunneling_directions.len());
        for (theta, phi, v) in tunneling_directions {
            println!("Direction (theta={:.4}, phi={:.4}) satisfies condition with f={:.6e}", theta, phi, v);
        }
    }
}
