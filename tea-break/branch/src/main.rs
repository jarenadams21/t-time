use std::f64::consts::PI;
use std::io::Write;

// Physical constants (CGS + suitable units):
const C: f64 = 3.0e10;        // speed of light in cm/s
const MEV_TO_ERG: f64 = 1.6021765e-6;
const KEV_TO_MEV: f64 = 1e-3;
const M_EC2: f64 = 0.511;     // electron rest energy in MeV
const R_E: f64 = 2.8179403262e-13; // classical electron radius (cm)

// Approximate scenario parameters:
const E_ISO_ERG: f64 = 1.0e55;      // isotropic energy in erg
const DURATION: f64 = 600.0;        // seconds
const R0: f64 = 1.0e13;             // initial radius in cm
const BETA: f64 = 1.0;              // ultra-relativistic approximation

// Band function parameters:
const ALPHA: f64 = -1.0;
const BETA_PAR: f64 = -2.3;
const E0_KEV: f64 = 300.0; // break energy in keV

struct Params {
    norm: f64, // normalization for the Band function
}

// Implement Clone for State so we can use `.clone()`
#[derive(Clone)]
struct State {
    time: f64,
    radius: f64,
    n_pairs: f64,
    n_photon: f64,
}

// Band function:
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
    let guess = 1.0;
    let mut params = Params{norm: guess};
    let total = total_energy_band(&params);
    let target = E_ISO_ERG;
    let adjusted = guess*(target/total);
    adjusted
}

fn pair_xsec_approx(E1_mev: f64, E2_mev: f64) -> f64 {
    let threshold = 4.0*M_EC2*M_EC2;
    let product = E1_mev*E2_mev;
    if product > threshold {
        PI*(R_E*R_E)
    } else {
        0.0
    }
}

fn pair_production_rate(n_ph: f64, params: &Params) -> f64 {
    let e_min = 1.0;
    let e_max = 1e5;
    let steps = 50;
    let de = (e_max - e_min)/(steps as f64);

    let mut spectrum = Vec::with_capacity(steps);
    let mut n_total = 0.0;

    for i in 0..steps {
        let E_kev = e_min + (i as f64)*de;
        let val = band_spectrum(E_kev, ALPHA, BETA_PAR, E0_KEV, params.norm);
        spectrum.push((E_kev, val));
        n_total += val*de;
    }

    for i in 0..steps {
        spectrum[i].1 /= n_total; // normalize to 1
    }

    let mut pair_rate = 0.0;
    for i in 0..steps {
        for j in 0..steps {
            let (E1_kev, f1) = spectrum[i];
            let (E2_kev, f2) = spectrum[j];
            let E1_mev = E1_kev*KEV_TO_MEV;
            let E2_mev = E2_kev*KEV_TO_MEV;
            let sigma = pair_xsec_approx(E1_mev, E2_mev);
            // Very rough dimension treatment:
            pair_rate += f1*f2*sigma*C*de*de;
        }
    }

    // scale by n_ph^2
    pair_rate*n_ph*n_ph
}

// Derivatives function:
fn derivatives(s: &State, params: &Params) -> f64 {
    pair_production_rate(s.n_photon, params)
}

// RK4 integrator:
fn rk4_step(s: &mut State, dt: f64, params: &Params) {
    let s_original = s.clone();
    let k1 = derivatives(&s_original, params);

    let mut s2 = s_original.clone();
    s2.time = s_original.time + dt/2.0;
    s2.n_pairs = s_original.n_pairs + k1*(dt/2.0);

    let k2 = derivatives(&s2, params);

    let mut s3 = s_original.clone();
    s3.time = s_original.time + dt/2.0;
    s3.n_pairs = s_original.n_pairs + k2*(dt/2.0);

    let k3 = derivatives(&s3, params);

    let mut s4 = s_original.clone();
    s4.time = s_original.time + dt;
    s4.n_pairs = s_original.n_pairs + k3*dt;

    let k4 = derivatives(&s4, params);

    s.n_pairs += (k1 + 2.0*k2 + 2.0*k3 + k4)*(dt/6.0);
}

impl State {
    fn new(n_photon_init: f64) -> Self {
        State {
            time: 0.0,
            radius: R0,
            n_pairs: 0.0,
            n_photon: n_photon_init,
        }
    }
}

fn main() {
    let norm = find_norm_for_band();
    let params = Params {norm};
    println!("Normalization for Band function: {}", params.norm);

    let avg_E_kev = 300.0;
    let avg_E_mev = avg_E_kev*KEV_TO_MEV;
    let avg_E_erg = avg_E_mev*MEV_TO_ERG;
    let vol = (4.0/3.0)*PI*R0.powi(3);
    let total_photons = E_ISO_ERG / avg_E_erg;
    let n_photon_init = total_photons/vol;

    let mut state = State::new(n_photon_init);

    let total_time = 3600.0;
    let dt = 1.0;

    let mut file = std::fs::File::create("simulation_output.csv").unwrap();
    writeln!(file, "time(s),radius(cm),n_photon(cm^-3),n_pairs(cm^-3)").unwrap();

    for _ in 0..(total_time as usize) {
        state.time += dt;
        state.radius = R0 + BETA*C*state.time;
        let scale = (R0/state.radius).powi(3);
        state.n_photon = n_photon_init*scale;

        rk4_step(&mut state, dt, &params);

        writeln!(file, "{},{},{},{}",
                 state.time,
                 state.radius,
                 state.n_photon,
                 state.n_pairs).unwrap();
    }

    println!("Final pairs: {} cm^-3", state.n_pairs);
    println!("Data in simulation_output.csv");
}