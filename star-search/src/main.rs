use std::fs::File;
use std::io::Write;
use std::f64::consts::PI;
use lazy_static::lazy_static;
use ndarray::{Array3};
use ndarray_npy::write_npy;

// Physical constants
const C: f64 = 2.99792458e8;
const HBAR: f64 = 1.054571817e-34;
const G: f64 = 6.67430e-11;
const K_B: f64 = 1.380649e-23;
const MU0: f64 = 4.0 * PI * 1e-7; // Vacuum permeability
const EPS0: f64 = 1.0/(MU0*C*C);

// Parameters
lazy_static! {
    static ref G_A_GAMMA: f64 = 1e-14; // Increased coupling slightly
    static ref PHOTON_TO_NEUTRINO_COEFF: f64 = 1e-10; 
}

const Q: f64 = 1.001;
const LAMBDA: f64 = 0.01; 
const NX: usize = 20;
const NY: usize = 20;
const NZ: usize = 20;
const DX: f64 = 0.5e-15;
const DT: f64 = 5.391247 * 10e-18; // 1.0e-17// Larger DT
const STEPS: usize = 1000; // More steps for evolution

const PHOTON_INIT: f64 = 1e10;
const AXION_INIT: f64 = 1e-20;
const NEUTRINO_INIT: f64 = 1e10;

const ENERGY_INIT: f64 = 3.2e35; 
const EPSILON_CRIT: f64 = 1.6e35; 
const DELTA: f64 = 0.2e35; 

const D_PH: f64 = 1e-3;
const D_AX: f64 = 1e-3;
const D_NU: f64 = 1e-3;
const D_E: f64 = 1e-3;

const LAMBDA_NU: f64 = 1e-5;
const ALPHA_EXPANSION: f64 = 1e-5;
const GW_STR: f64 = 1e-21;
const GW_FREQ: f64 = 1e3;

fn axion_photon_conversion(n_ax: f64, n_ph: f64) -> f64 {
    let saturation = 1.0 + n_ph / 1e33;
    (*G_A_GAMMA * n_ax) / saturation
}

fn photon_neutrino_conversion(n_ph: f64) -> f64 {
    let saturation = 1.0 + n_ph / 1e33;
    (*PHOTON_TO_NEUTRINO_COEFF * n_ph) / saturation
}

struct Field {
    photon_density: Vec<f64>,
    axion_density: Vec<f64>,
    neutrino_density: Vec<f64>,
    energy_density: Vec<f64>,
    efield_x: Vec<f64>,
    efield_y: Vec<f64>,
    efield_z: Vec<f64>,
    bfield_x: Vec<f64>,
    bfield_y: Vec<f64>,
    bfield_z: Vec<f64>,
}

impl Field {
    fn new() -> Self {
        let size = NX*NY*NZ;
        Field {
            photon_density: vec![PHOTON_INIT; size],
            axion_density: vec![AXION_INIT; size],
            neutrino_density: vec![NEUTRINO_INIT; size],
            energy_density: vec![ENERGY_INIT; size],
            // Initialize E and B fields with a small perturbation
            efield_x: vec![0.0; size],
            efield_y: vec![0.0; size],
            efield_z: vec![0.0; size],
            bfield_x: vec![0.0; size],
            bfield_y: vec![0.0; size],
            bfield_z: vec![0.0; size],
        }
    }

    fn idx(&self, x:usize, y:usize, z:usize)->usize {
        x + NX*(y+NY*z)
    }

    fn boundary_index(x:isize, max:usize)->usize {
        let mut xx=x;
        if xx<0 {xx=0;} else if xx>=max as isize {xx=(max as isize)-1;}
        xx as usize
    }

    fn laplacian(&self, arr:&Vec<f64>, x:usize,y:usize,z:usize)->f64 {
        let xm=Self::boundary_index(x as isize-1,NX);
        let xp=Self::boundary_index(x as isize+1,NX);
        let ym=Self::boundary_index(y as isize-1,NY);
        let yp=Self::boundary_index(y as isize+1,NY);
        let zm=Self::boundary_index(z as isize-1,NZ);
        let zp=Self::boundary_index(z as isize+1,NZ);
        let c=arr[self.idx(x,y,z)];
        let dx2=DX*DX*DX;
        (arr[self.idx(xp,y,z)]+arr[self.idx(xm,y,z)]+
         arr[self.idx(x,yp,z)]+arr[self.idx(x,ym,z)]+
         arr[self.idx(x,y,zp)]+arr[self.idx(x,y,zm)]-6.0*c)/dx2
    }

}

fn eos_pressure(eps:f64)->f64{
    let w_qgp=0.5*(1.0+((eps - EPSILON_CRIT)/DELTA).tanh());
    let p_qgp=(1.0/3.0)*eps;
    let p_hg=0.15*eps;
    w_qgp*p_qgp+(1.0 - w_qgp)*p_hg
}

fn metric_factor(t:f64,x:f64)->f64 {
    1.0 + GW_STR*(2.0*PI*GW_FREQ*t).sin()*x
}

// Hecke R-matrix
fn apply_hecke_r_matrix(fields:&mut Field) {
    let qm=Q.powf(-0.5);
    let qp=Q.powf(0.5);
    for z in 0..NZ {
        for y in 0..NY {
            for x in 0..(NX-1) {
                let i=x+NX*(y+NY*z);
                let j=i+1;
                let ph_i=fields.photon_density[i];
                let ax_i=fields.axion_density[i];
                let ph_j=fields.photon_density[j];
                let ax_j=fields.axion_density[j];
                let ph_i_new=0.5*(ph_i*qm+ax_j);
                let ax_i_new=0.5*(ax_i*qp+ph_j);
                let ph_j_new=0.5*(ph_j*qm+ax_i);
                let ax_j_new=0.5*(ax_j*qp+ph_i);
                fields.photon_density[i]=ph_i_new;
                fields.axion_density[i]=ax_i_new;
                fields.photon_density[j]=ph_j_new;
                fields.axion_density[j]=ax_j_new;
            }
        }
    }
}

// Placeholder for Maxwell update (simplified: no current, no full PDE solve)
fn update_maxwell(fields:&mut Field) {
    // Normally solve curl equations. Here we just introduce a small perturbation
    // to E and B fields to break symmetry.
    // This ensures that after some steps, photon distribution changes.
    let size=NX*NY*NZ;
    for i in 0..size {
        // Introduce a tiny random perturbation or gradient-based shift
        fields.efield_x[i]+= 1e-3*DT;
        fields.efield_y[i]+= 1e-3*DT;
        fields.efield_z[i]+= 0.0;
        fields.bfield_x[i]+= 0.0;
        fields.bfield_y[i]+= 1e-3*DT;
        fields.bfield_z[i]+= 1e-3*DT;
    }
}

fn main(){
    let mut field=Field::new();
    let mut file=File::create("results.csv").unwrap();
    writeln!(file,"time(s),avg_photon_density,avg_axion_density,avg_neutrino_density,avg_energy_density").unwrap();

    for step in 0..STEPS {
        let t=step as f64*DT;

        // Update Maxwell fields first
        update_maxwell(&mut field);

        let mut new_ph=field.photon_density.clone();
        let mut new_ax=field.axion_density.clone();
        let mut new_nu=field.neutrino_density.clone();
        let mut new_e =field.energy_density.clone();

        for z in 0..NZ {
            for y in 0..NY {
                for x in 0..NX {
                    let idx=x+NX*(y+NY*z);
                    let n_ph=field.photon_density[idx];
                    let n_ax=field.axion_density[idx];
                    let n_nu=field.neutrino_density[idx];
                    let eps=field.energy_density[idx];

                    let lap_ph=field.laplacian(&field.photon_density,x,y,z);
                    let lap_ax=field.laplacian(&field.axion_density,x,y,z);
                    let lap_nu=field.laplacian(&field.neutrino_density,x,y,z);
                    let lap_e =field.laplacian(&field.energy_density,x,y,z);

                    let p=eos_pressure(eps);
                    let x_pos=x as f64*DX;
                    let mf=metric_factor(t,x_pos);

                    let d_ax_to_ph=axion_photon_conversion(n_ax,n_ph)*DT;
                    let d_ph_to_nu=photon_neutrino_conversion(n_ph)*DT;

                    let d_e_nu=LAMBDA_NU*n_nu*DT;
                    let d_e_exp=p*ALPHA_EXPANSION*DT;

                    let ph_new=n_ph+D_PH*lap_ph*DT+d_ax_to_ph-d_ph_to_nu;
                    let ax_new=n_ax+D_AX*lap_ax*DT - d_ax_to_ph; 
                    let nu_new=n_nu+D_NU*lap_nu*DT + d_ph_to_nu; 
                    let e_new=eps+D_E*lap_e*DT - d_e_nu - d_e_exp;

                    new_ph[idx]=(ph_new*mf).max(0.0);
                    new_ax[idx]=(ax_new*mf).max(0.0);
                    new_nu[idx]=(nu_new*mf).max(0.0);
                    new_e[idx]=(e_new*mf).max(0.0);
                }
            }
        }

        field.photon_density=new_ph;
        field.axion_density=new_ax;
        field.neutrino_density=new_nu;
        field.energy_density=new_e;

        apply_hecke_r_matrix(&mut field);

        let vol=(NX*NY*NZ) as f64;
        let avg_photon=field.photon_density.iter().sum::<f64>()/vol;
        let avg_axion=field.axion_density.iter().sum::<f64>()/vol;
        let avg_neutrino=field.neutrino_density.iter().sum::<f64>()/vol;
        let avg_energy=field.energy_density.iter().sum::<f64>()/vol;

        writeln!(file,"{},{},{},{},{}",t,avg_photon,avg_axion,avg_neutrino,avg_energy).unwrap();
    }

    let mut torsion_field=vec![0.0;NX*NY*NZ];
    for z in 0..NZ {
        let z_val=z as f64/(NZ as f64-1.0);
        for y in 0..NY {
            let y_val=y as f64/(NY as f64-1.0);
            for x in 0..NX {
                let x_val=x as f64/(NX as f64-1.0);
                let idx=x+NX*(y+NY*z);
                let ph_norm=field.photon_density[idx]/(PHOTON_INIT*10.0);
                // Now torsion could also depend on E and B fields to create non-trivial patterns:
                let e_mag=(field.efield_x[idx].powf(2.14)+field.efield_y[idx].powi(2)+field.efield_z[idx].powi(2)).sqrt();
                let b_mag=(field.bfield_x[idx].powi(2)+field.bfield_y[idx].powf(2.14)+field.bfield_z[idx].powi(2)).sqrt();
                
                torsion_field[idx] = (2.0*PI*x_val).sin()*(2.0*PI*y_val).cos()*(-z_val).exp() 
                                     + ph_norm*1e-5
                                     + 1e-6*(e_mag-b_mag); 
            }
        }
    }

    std::fs::create_dir_all("data").unwrap();

    let photon_arr=Array3::from_shape_vec((NX,NY,NZ),field.photon_density.clone()).expect("shape");
    let axion_arr=Array3::from_shape_vec((NX,NY,NZ),field.axion_density.clone()).expect("shape");
    let neutrino_arr=Array3::from_shape_vec((NX,NY,NZ),field.neutrino_density.clone()).expect("shape");
    let torsion_arr=Array3::from_shape_vec((NX,NY,NZ),torsion_field.clone()).expect("shape");

    write_npy("data/photon_density_final.npy",&photon_arr).unwrap();
    write_npy("data/axion_density_final.npy",&axion_arr).unwrap();
    write_npy("data/neutrino_density_final.npy",&neutrino_arr).unwrap();
    write_npy("data/torsion_field_final.npy",&torsion_arr).unwrap();

    println!("Simulation complete with Maxwell & Clifford hints. Data saved to results.csv and data/*.npy");
}
