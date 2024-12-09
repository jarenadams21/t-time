use std::fs::File;
use std::io::Write;

// We assume a monoatomic ideal gas with degrees of freedom f=3, γ = Cp/Cv = 5/3 for demonstration.
// Dimensionless constants: k_B = 1, N = 1, so P V = T applies.

// Carnot cycle parameters
const T_H: f64 = 425.0;   // Hot reservoir temperature (K)
const T_C: f64 = 278.0;   // Cold reservoir temperature (K)
const GAMMA: f64 = 5.0/3.0; // Ratio of specific heats for monoatomic ideal gas
const N_STEPS_ISOTHERMAL: usize = 500; // Number of steps for isothermal processes
const N_STEPS_ADIABATIC: usize = 500;  // Number of steps for adiabatic processes

// Initial state A: T_H and choose V1=1.0
// From ideal gas: P1 = T_H / V1 = T_H since V1=1.0
const V1: f64 = 1.0;
const T_START: f64 = T_H;

fn main() {
    let mut file = File::create("carnot_data.csv").unwrap();
    writeln!(file, "step,V,P,T,S,Q,W,phase").unwrap();

    // We pick a convenient ratio for the first isothermal expansion: V2 = 2.0 * V1
    let v2 = 2.0;

    // From B to C (adiabatic): T_H * V2^(γ-1) = T_C * V3^(γ-1)
    // => V3 = V2 * (T_H/T_C)^(1/(γ-1))
    let ratio = (T_H/T_C).powf(1.0/(GAMMA-1.0));
    let v3 = v2 * ratio;

    // From C to D (isothermal at T_C), we compress down to V4
    // From D to A (adiabatic): T_C * V4^(γ-1) = T_H * V1^(γ-1)
    // => V4 = (T_H/T_C)^(1/(γ-1)) ≈ 2.15
    let v4 = ratio;

    // The cycle is: A(500K,1.0), B(500K,2.0), C(300K,4.3), D(300K,2.15), A(500K,1.0)
    // Perfectly closes the loop.

    // Initialize
    let mut v = V1;
    let mut t = T_START;
    let s_ref = v.ln() + 1.5 * t.ln(); // Entropy reference at state A
    let mut q_cumulative = 0.0;
    let mut w_cumulative = 0.0;
    let mut step_count = 0;

    let mut write_data = |step: usize, v_val: f64, t_val: f64, q_val: f64, w_val: f64, phase: &mut File, s_ref: f64| {
        let p_val = t_val / v_val;
        let s_val = v_val.ln() + 1.5 * t_val.ln() - s_ref;
        writeln!(phase, "{},{},{},{},{},{},{},{}",
                 step, v_val, p_val, t_val, s_val, q_val, w_val, "phase").unwrap();
    };

    // Helper with explicit phase name
    let mut write_data_phase = |step: usize, v_val: f64, t_val: f64, q_val: f64, w_val: f64, phase_name: &str, file: &mut File, s_ref: f64| {
        let p_val = t_val / v_val;
        let s_val = v_val.ln() + 1.5 * t_val.ln() - s_ref;
        writeln!(file, "{},{},{},{},{},{},{},{}",
                 step, v_val, p_val, t_val, s_val, q_val, w_val, phase_name).unwrap();
    };

    // 1) Isothermal expansion at T_H: A->B
    // V: 1.0 to 2.0 in N_STEPS_ISOTHERMAL
    let dv_iso_hot = (v2 - v) / (N_STEPS_ISOTHERMAL as f64);
    for _ in 0..N_STEPS_ISOTHERMAL {
        let p_local = t / v;
        let dW = p_local * dv_iso_hot;
        let dQ = dW; // isothermal => ΔU=0 => Q=W
        w_cumulative += dW;
        q_cumulative += dQ;
        v += dv_iso_hot;
        step_count += 1;
        write_data_phase(step_count, v, t, q_cumulative, w_cumulative, "isothermal_hot", &mut file, s_ref);
    }

    // 2) Adiabatic expansion B->C: V: 2.0 to V3
    // Adiabatic: T * V^(γ-1) = const
    let c_adiab1 = t * v.powf(GAMMA - 1.0);
    let dv_adiab_expand = (v3 - v) / (N_STEPS_ADIABATIC as f64);
    for _ in 0..N_STEPS_ADIABATIC {
        v += dv_adiab_expand;
        t = c_adiab1 / v.powf(GAMMA - 1.0);
        let p_local = t / v;
        let dW = p_local * dv_adiab_expand;
        w_cumulative += dW;
        // Q=0 adiabatic
        step_count += 1;
        write_data_phase(step_count, v, t, q_cumulative, w_cumulative, "adiabatic_expand", &mut file, s_ref);
    }

    // 3) Isothermal compression at T_C: C->D
    // V: v3 to v4 at T_C
    t = T_C;
    let mut v_current = v3;
    let dv_iso_cold = (v4 - v3) / (N_STEPS_ISOTHERMAL as f64);
    for _ in 0..N_STEPS_ISOTHERMAL {
        let p_local = t / v_current;
        let dW = p_local * dv_iso_cold;
        // Compression: from system perspective, dW<0. Q=W for isothermal
        w_cumulative += dW;
        q_cumulative += dW;
        v_current += dv_iso_cold;
        step_count += 1;
        write_data_phase(step_count, v_current, t, q_cumulative, w_cumulative, "isothermal_cold", &mut file, s_ref);
    }

    v = v_current;

    // 4) Adiabatic compression D->A: V: v4 to V1 at final T_H
    let c_adiab2 = t * v.powf(GAMMA - 1.0);
    let dv_adiab_back = (V1 - v) / (N_STEPS_ADIABATIC as f64);
    for _ in 0..N_STEPS_ADIABATIC {
        v += dv_adiab_back;
        t = c_adiab2 / v.powf(GAMMA - 1.0);
        let p_local = t / v;
        let dW = p_local * dv_adiab_back;
        w_cumulative += dW;
        // Q=0 adiabatic
        step_count += 1;
        write_data_phase(step_count, v, t, q_cumulative, w_cumulative, "adiabatic_compress", &mut file, s_ref);
    }

    println!("Simulation complete. Data in carnot_data.csv");
}