use std::f64::consts::PI;
use std::io::Write;

// Physical constants
const C: f64 = 2.99792458e8;       // m/s
const HBAR: f64 = 1.054571817e-34; // J·s
const K_B: f64 = 1.380649e-23;     // J/K
const G: f64 = 6.67430e-11;        // m^3 kg^-1 s^-2
const H_0: f64 = 67.4e3/(3.086e22);// Hubble constant ~ 2.2e-18 s^-1 (67.4 km/s/Mpc)
const SIGMA_T: f64 = 6.6524587321e-29; // Thomson cross section

// Cosmological parameters today (for simplicity):
const OMEGA_M: f64 = 0.315;
const OMEGA_R: f64 = 9.0e-5; // radiation today
const OMEGA_L: f64 = 1.0 - OMEGA_M - OMEGA_R; // flat Universe

// CMB temperature today:
const T0: f64 = 2.7255; // K
// Photon number density today ~ 4.11e8 m^-3
const N_GAMMA_0: f64 = 4.11e8; 
// Assume a primordial B-field upper limit:
const B0: f64 = 1e-9; // Tesla
// Axion-photon coupling upper limit:
const G_AGAMMA: f64 = 1e-20; // J^-1 approx
// Choose a box size and grid:
const NX: usize = 20;
const NY: usize = 20;
const NZ: usize = 20;
const DX: f64 = 1.0; // m
const DT: f64 = 1e11; // s (large time step to simulate cosmic evolution)
const TOTAL_TIME: f64 = 4.35e17; // ~13.8 Gyr in seconds

struct Field3D {
    photons: Vec<f64>,
    axions: Vec<f64>,
}

impl Field3D {
    fn new(ph_init: f64) -> Self {
        let size = NX*NY*NZ;
        Field3D {
            photons: vec![ph_init; size],
            axions: vec![0.0; size],
        }
    }

    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + NX * (y + NY * z)
    }

    fn laplacian(&self, arr: &Vec<f64>, x: usize, y: usize, z: usize) -> f64 {
        let val = arr[self.idx(x,y,z)];
        let xp = arr[self.idx((x+1)%NX,y,z)];
        let xm = arr[self.idx((x+NX-1)%NX,y,z)];
        let yp = arr[self.idx(x,(y+1)%NY,z)];
        let ym = arr[self.idx(x,(y+NY-1)%NY,z)];
        let zp = arr[self.idx(x,y,(z+1)%NZ)];
        let zm = arr[self.idx(x,y,(z+NZ-1)%NZ)];
        (xp+xm+yp+ym+zp+zm - 6.0*val)/(DX*DX)
    }
}

// Friedmann equation solver for a(t):
// da/dt = a * H(a), with H(a) from LCDM:
fn hubble(a: f64) -> f64 {
    H_0 * (OMEGA_R/a.powi(4) + OMEGA_M/a.powi(3) + OMEGA_L).sqrt()
}

// Axion conversion rate:
fn axion_rate(a: f64) -> f64 {
    // B(t) = B0*(a0/a)^2 with a0=1 today
    let b = B0/(a*a);
    // On scale L = DX, P ~ (g_{aγ} B L / (ħc))^2 per segment of coherence
    let p = (G_AGAMMA*b*DX/(HBAR*C)).powi(2);
    // rate ~ p*c/L to get transitions per second:
    (p*C/DX)
}

// Photon diffusion coefficient (scaled):
// n_e ~ depends on reionization, etc., assume n_e negligible today?
// For demonstration, take n_e ~ small. Realistically, after reionization:
// n_e ~ 0.5 * n_baryon ~ 0.5 * (current baryon density ~0.22/m^3) ~0.1/m^3 (very rough)
// This is extremely small, meaning scattering negligible today.
// We'll set D=0 for modern era, just to show the form:
fn diffusion_coefficient(_a: f64) -> f64 {
    // For realistic present universe, photon free path is huge, D large.
    // On small scale, negligible. Just return 0 here.
    0.0
}

fn main() {
    // initial conditions: start from early universe: set a start at a ~ 1e-3 (z~999)
    let mut a = 1e-3; 
    let mut t = 0.0;

    // Photon number density at scale factor a: n_gamma = N_GAMMA_0 / a^3
    let n_ph_init = N_GAMMA_0/(a*a*a);
    let mut field = Field3D::new(n_ph_init);

    let mut file = std::fs::File::create("cosmic_evolution.csv").unwrap();
    writeln!(file, "time(s),scale_factor,a,avg_photon(m^-3),avg_axion(m^-3)").unwrap();

    while t < TOTAL_TIME {
        // Compute Hubble rate and evolve a(t):
        let H = hubble(a);
        // da/dt = a * H
        let da = a * H * DT;
        a += da;
        t += DT;

        // Update fields:
        let ax_rate = axion_rate(a);
        let D = diffusion_coefficient(a);

        // Evolve photon and axion fields:
        let mut new_photons = field.photons.clone();
        let mut new_axions = field.axions.clone();
        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx = field.idx(x,y,z);
                    let n_ph = field.photons[idx];
                    let n_ax = field.axions[idx];

                    // Redshift scaling (if we were to do small steps over cosmic time, we'd adjust density)
                    // Over a short DT, the change in n_ph from expansion: n_ph ~ n_ph / a^3, 
                    // Let's apply n_ph_new = n_ph*(a_old^3 / a_new^3):
                    // But we are evolving a in large steps. For stability, do tiny increments or store old a:
                    // Approximate scaling from last step:
                    let a_old = a - da;
                    let scale_factor_ratio = (a_old/a).powi(3);
                    let n_ph_expanded = n_ph * scale_factor_ratio;

                    // Diffusion (small scale - likely negligible now)
                    let lap = field.laplacian(&field.photons,x,y,z);
                    let dn_ph_diff = D * lap * DT;

                    // Axion production:
                    // dn_ax ~ n_ph * axion_rate * DT
                    let dn_ax = n_ph_expanded * ax_rate * DT;

                    let new_n_ph = n_ph_expanded + dn_ph_diff - dn_ax;
                    let new_n_ax = n_ax * scale_factor_ratio + dn_ax; // Axions also diluted by expansion

                    new_photons[idx] = new_n_ph;
                    new_axions[idx] = new_n_ax;
                }
            }
        }

        field.photons = new_photons;
        field.axions = new_axions;

        let avg_ph = field.photons.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        let avg_ax = field.axions.iter().sum::<f64>()/(NX*NY*NZ) as f64;
        writeln!(file, "{},{},{},{}", t, a, avg_ph, avg_ax).unwrap();
    }

    println!("Simulation completed. Results in cosmic_evolution.csv");
}
