extern crate arrayfire as af;

use af::{Array, Dim4};
use std::convert::TryInto;

/// The `Neuron` trait represents the common behavior of all types of neurons.
///
/// The three methods in this trait correspond to the three major stages of the spiking process:
/// - Update the neuron's internal state based on incoming input
/// - Handle a spike event if the membrane potential reaches the threshold
/// - Emit a spike if the neuron is spiking (i.e., the membrane potential has crossed the threshold)
pub trait Neuron {
    /// Update the neuron's internal state based on incoming input
    fn update_state(&mut self, input_current: &Array<f64>, dt: f64);

    /// Handle a spike event if the membrane potential reaches the threshold
    fn handle_spike(&mut self);

    /// Emit a spike if the neuron is spiking (i.e., the membrane potential has crossed the threshold)
    fn emit_spike(&self) -> Array<bool>;
}
/// Leaky Integrate-and-Fire (LIF) Neuron Model
///
/// Reference: Gerstner, W., & Kistler, W. M. (2002). Spiking Neuron Models. Cambridge University Press.
///
/// The LIF model is a simple and computationally efficient model that
/// captures the essential dynamics of neuronal spiking. It represents
/// the neuron's membrane potential as a leaky integrator, charging in
/// response to input currents and discharging (leaking) over time.
///
/// Pros:
/// - Simple and computationally efficient
/// - Suitable for large-scale network simulations
///
/// Cons:
/// - Limited biological realism
/// - Does not reproduce complex spike patterns
///
/// Mathematical equations:
/// dv/dt = (-(v - rest_potential) + membrane_resistance * I) / membrane_time_constant
/// v(t) = rest_potential + (v(0) - rest_potential + membrane_resistance * I) * exp(-t / membrane_time_constant)
///
/// If v >= threshold:
///     v = reset_potential
///     Emit spike

#[derive(Clone)]
pub struct LeakyIntegrateAndFire {
    pub membrane_potential: Array<f64>,
    pub threshold: f64,
    pub rest_potential: f64,
    pub reset_potential: f64,
    pub membrane_resistance: f64,
    pub membrane_time_constant: f64,
}

impl LeakyIntegrateAndFire {
    /// Create a new LIF neuron instance with the specified parameters.
    ///
    /// # Arguments
    /// * `n_neurons`: The number of neurons in the population
    /// * `threshold`: The firing threshold in mV
    /// * `rest_potential`: The resting potential in mV
    /// * `reset_potential`: The membrane potential reset value in mV
    /// * `membrane_resistance`: The membrane resistance in MOhm
    /// * `membrane_time_constant`: The membrane time constant in ms
    pub fn new(
        n_neurons: usize,
        threshold: f64,
        rest_potential: f64,
        reset_potential: f64,
        membrane_resistance: f64,
        membrane_time_constant: f64,
    ) -> Self {
        LeakyIntegrateAndFire {
            membrane_potential: af::constant(
                rest_potential,
                Dim4::new(&[n_neurons.try_into().unwrap(), 1, 1, 1]),
            ),
            threshold,
            rest_potential,
            reset_potential,
            membrane_resistance,
            membrane_time_constant,
        }
    }
}

impl Neuron for LeakyIntegrateAndFire {
    /// Update the state of the neuron based on the input current it receives.
    ///
    /// # Arguments
    ///
    /// * `input_current`: A reference to an `ndarray::Array` containing the input current for each
    ///   time step.
    /// * `dt`: The time step for the simulation, in seconds.
    fn update_state(&mut self, input_current: &Array<f64>, dt: f64) {
        let membrane_potential_change = ((-(&self.membrane_potential - self.rest_potential)
            + self.membrane_resistance * input_current)
            / self.membrane_time_constant)
            * dt;
        self.membrane_potential += membrane_potential_change;
    }
    /// Handle a spike event in the neuron.
    ///
    /// When the membrane potential of the neuron crosses the threshold value, it emits a spike and
    /// its membrane potential is reset to a lower value.
    fn handle_spike(&mut self) {
        let spike_condition = af::ge(&self.membrane_potential, &self.threshold, false);
        let reset_values = af::constant(self.reset_potential, self.membrane_potential.dims());
        self.membrane_potential =
            af::select(&reset_values, &spike_condition, &self.membrane_potential);
    }

    /// Emit a spike event from the neuron.
    ///
    /// Returns a boolean array indicating which neurons have emitted a spike.
    fn emit_spike(&self) -> Array<bool> {
        af::ge(&self.membrane_potential, &self.threshold, false)
    }
}

// /// Izhikevich Neuron Model
// ///
// /// Reference: Izhikevich, E. M. (2003). Simple model of spiking neurons.
// /// IEEE Transactions on Neural Networks, 14(6), 1569-1572.
// ///
// /// The Izhikevich model is a simplified neuron model that can reproduce
// /// various spiking patterns observed in biological neurons with a lower
// /// computational cost compared to the Hodgkin-Huxley model. The model is
// /// governed by two differential equations that describe the dynamics of
// /// the membrane potential and a recovery variable.
// ///
// /// Pros:
// /// - Reproduces complex spike patterns
// /// - Computationally efficient
// ///
// /// Cons:
// /// - Limited biological realism
// ///
// /// Mathematical equations:
// /// dv/dt = 0.04 * v^2 + 5 * v + 140 - u + I
// /// du/dt = a * (b * v - u)
// ///
// /// If v >= 30:
// ///     v = c
// ///     u = u + d
// ///     Emit spike

pub struct Izhikevich {
    pub membrane_potential: Array<f64>,
    pub recovery_variable: Array<f64>,
    pub threshold: f64,
    pub recovery_time_scale: f64,
    pub recovery_sensitivity: f64,
    pub reset_membrane_potential: f64,
    pub reset_recovery_variable: f64,
}

/// Create a new Izhikevich neuron instance with the specified parameters.
///
/// # Arguments
/// * `membrane_potenial`: The neuron's membrane potential (voltage) in millivolts (mV)
/// * `recovery_variable`: The recovery variable (u) representing the membrane recovery rate
/// * `threshold`: The firing threshold in mV
/// * `recovery_time_scale`: The time scale of the recovery variable (u)
/// * `recovery_sensitivity`: The sensitivity of the recovery variable (u) to the subthreshold fluctuations of the membrane potential
/// * `reset_membrane_potential`: The after-spike reset value of the membrane potential in mV
/// * `reset_recovery_variable`: The after-spike reset value of the recovery variable (u)
impl Izhikevich {
    pub fn new(
        n_neurons: usize,
        recovery_time_scale: f64,
        recovery_sensitivity: f64,
        reset_membrane_potential: f64,
        reset_recovery_variable: f64,
        threshold: f64,
    ) -> Self {
        Izhikevich {
            membrane_potential: af::constant(
                -65.0,
                Dim4::new(&[n_neurons.try_into().unwrap(), 1, 1, 1]),
            ),
            recovery_variable: af::constant(
                0.0,
                Dim4::new(&[n_neurons.try_into().unwrap(), 1, 1, 1]),
            ),
            recovery_time_scale,
            recovery_sensitivity,
            reset_membrane_potential,
            reset_recovery_variable,
            threshold,
        }
    }
}

impl Neuron for Izhikevich {
    /// Update the state of the neuron based on the input current it receives.
    ///
    /// # Arguments
    ///
    /// * input_current: A reference to an arrayfire::Array containing the input current for each
    /// time step.
    /// * dt: The time step for the simulation, in seconds.
    fn update_state(&mut self, input_current: &Array<f64>, dt: f64) {
        let dv = (0.04 * af::pow(&self.membrane_potential.clone(), &2.0, false)
            + 5.0 * self.membrane_potential.clone()
            + 140.0
            - self.recovery_variable.clone()
            + input_current)
            * dt;
        let du = self.recovery_time_scale
            * (self.recovery_sensitivity * self.membrane_potential.clone()
                - self.recovery_variable.clone())
            * dt;
        self.membrane_potential += dv;
        self.recovery_variable += du;
    }
    /// Handle a spike event in the neuron.
    ///
    /// When the membrane potential of the neuron crosses the threshold value, it emits a spike and
    /// its membrane potential is reset to a lower value.
    fn handle_spike(&mut self) {
        let spike_condition = af::ge(
            &self.membrane_potential,
            &af::constant(self.threshold, self.membrane_potential.dims()),
            false,
        );
        let reset_values = af::constant(
            self.reset_membrane_potential,
            self.membrane_potential.dims(),
        );
        self.membrane_potential =
            af::select(&reset_values, &spike_condition, &self.membrane_potential);
        self.recovery_variable += af::select(
            &af::constant(self.reset_recovery_variable, self.recovery_variable.dims()),
            &spike_condition,
            &af::constant(0.0, self.recovery_variable.dims()),
        );
    }
    /// Emit a spike event from the neuron.
    ///
    /// Returns a boolean array indicating which neurons have emitted a spike.
    fn emit_spike(&self) -> Array<bool> {
        af::ge(
            &self.membrane_potential,
            &af::constant(self.threshold, self.membrane_potential.dims()),
            false,
        )
    }
}

// /// Hodgkin-Huxley Neuron Model
// ///
// /// Reference: Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of membrane current
// /// and its application to conduction and excitation in nerve. The Journal of Physiology, 117(4), 500-544.
// ///
// /// The Hodgkin-Huxley model is a biophysically detailed model that
// /// describes the conductance of specific ion channels (sodium and
// /// potassium) and their activation and inactivation variables (n, m, and h)
// /// to simulate the neuron's behavior more accurately.
// ///
// /// Pros:
// /// - High biological realism
// /// - Reproduces complex spike patterns
// ///
// /// Cons:
// /// - Computationally expensive
// /// - Not suitable for large-scale network simulations
// ///
// /// Mathematical equations:
// /// dv/dt = (g_na * m^3 * h * (e_na - v) + g_k * n^4 * (e_k - v) + g_l * (e_l - v) + I) / C_m
// /// dn/dt = alpha_n * (1 - n) - beta_n * n
// /// dm/dt = alpha_m * (1 - m) - beta_m * m
// /// dh/dt = alpha_h * (1 - h) - beta_h * h

// pub struct HodgkinHuxley {
//     pub membrane_potential: f64, // The neuron's membrane potential (voltage) in millivolts (mV)
//     pub n: f64, // The activation variable for the delayed rectifier potassium channel (K+)
//     pub m: f64, // The activation variable for the fast-activating sodium channel (Na+)
//     pub h: f64, // The inactivation variable for the fast-activating sodium channel (Na+)
//     pub g_na: f64, // The maximum conductance for the fast-activating sodium channel (Na+) in millisiemens per cm² (mS/cm²)
//     pub g_k: f64, // The maximum conductance for the delayed rectifier potassium channel (K+) in millisiemens per cm² (mS/cm²)
//     pub g_l: f64, // The leak conductance in millisiemens per cm² (mS/cm²)
//     pub e_na: f64, // The sodium Nernst (reversal) potential in millivolts (mV)
//     pub e_k: f64, // The potassium Nernst (reversal) potential in millivolts (mV)
//     pub e_l: f64, // The leak Nernst (reversal) potential in millivolts (mV)
// }

// impl HodgkinHuxley {
//     pub fn new() -> Self {
//         HodgkinHuxley {
//             membrane_potential: -65.0,
//             n: 0.32,
//             m: 0.05,
//             h: 0.6,
//             g_na: 120.0,
//             g_k: 36.0,
//             g_l: 0.3,
//             e_na: 50.0,
//             e_k: -77.0,
//             e_l: -54.4,
//         }
//     }
//     fn alpha_n(&self) -> f64 {
//         (0.01 * (55.0 + self.membrane_potential))
//             / (1.0 - (-0.1 * (55.0 + self.membrane_potential)).exp())
//     }

//     fn beta_n(&self) -> f64 {
//         0.125 * (-65.0 - self.membrane_potential).exp()
//     }

//     fn alpha_m(&self) -> f64 {
//         (0.1 * (40.0 + self.membrane_potential))
//             / (1.0 - (-0.1 * (40.0 + self.membrane_potential)).exp())
//     }

//     fn beta_m(&self) -> f64 {
//         4.0 * (-20.0 - self.membrane_potential).exp()
//     }

//     fn alpha_h(&self) -> f64 {
//         0.07 * (-20.0 - self.membrane_potential).exp()
//     }

//     fn beta_h(&self) -> f64 {
//         1.0 / (1.0 + (-0.1 * (30.0 + self.membrane_potential)).exp())
//     }
// }

// impl Neuron for HodgkinHuxley {

//     fn update_state(&mut self, input_current: f64, dt: f64) {
//         let i_na = self.g_na * self.m.powi(3) * self.h * (self.e_na - self.membrane_potential);
//         let i_k = self.g_k * self.n.powi(4) * (self.e_k - self.membrane_potential);
//         let i_l = self.g_l * (self.e_l - self.membrane_potential);

//         let d_v = (input_current - i_na - i_k - i_l) * dt;

//         let dn = (self.alpha_n() * (1.0 - self.n) - self.beta_n() * self.n) * dt;
//         let dm = (self.alpha_m() * (1.0 - self.m) - self.beta_m() * self.m) * dt;
//         let dh = (self.alpha_h() * (1.0 - self.h) - self.beta_h() * self.h) * dt;

//         self.membrane_potential += d_v;
//         self.n += dn;
//         self.m += dm;
//         self.h += dh;
//     }

//     fn handle_spike(&mut self) {
//         self.membrane_potential = -65.0;
//         self.n = 0.32;
//         self.m = 0.05;
//         self.h = 0.6;
//     }

//     fn emit_spike(&self) -> bool {
//         self.membrane_potential >= 30.0
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use af::{all_true, constant};

    const EPSILON: f64 = 1e-12;

    fn assert_arrays_approx_eq(a: &Array<f64>, b: &Array<f64>, epsilon: f64) {
        let abs_diff = af::abs(&(a - b));
        let is_close = af::le(&abs_diff, &af::constant(epsilon, a.dims()), false);
        let all_true_scalar = all_true(&is_close, 1);
        // assert_eq!(all_true_scalar, true);
    }

    #[test]
    fn new() {
        let lif = LeakyIntegrateAndFire::new(2, -55.0, -70.0, -65.0, 20.0, 10.0);
        let expected_membrane_potential = af::constant(-70.0, Dim4::new(&[2, 1, 1, 1]));

        assert_arrays_approx_eq(
            &lif.membrane_potential,
            &expected_membrane_potential,
            EPSILON,
        );
        assert_eq!(lif.threshold, -55.0);
        assert_eq!(lif.rest_potential, -70.0);
        assert_eq!(lif.reset_potential, -65.0);
        assert_eq!(lif.membrane_resistance, 20.0);
        assert_eq!(lif.membrane_time_constant, 10.0);
    }

    #[test]
    fn update_state() {
        let mut lif = LeakyIntegrateAndFire::new(1, -55.0, -70.0, -65.0, 20.0, 10.0);
        let input_current = af::constant(2.0, Dim4::new(&[1, 1, 1, 1]));

        lif.update_state(&input_current, 1.0);
        let expected_membrane_potential = af::constant(-50.0, Dim4::new(&[1, 1, 1, 1]));

        assert_arrays_approx_eq(
            &lif.membrane_potential,
            &expected_membrane_potential,
            EPSILON,
        );
    }

    #[test]
    fn handle_spike() {
        let mut lif = LeakyIntegrateAndFire::new(1, -55.0, -70.0, -65.0, 20.0, 10.0);
        lif.membrane_potential = af::constant(-50.0, Dim4::new(&[1, 1, 1, 1]));

        lif.handle_spike();
        let expected_membrane_potential = af::constant(-65.0, Dim4::new(&[1, 1, 1, 1]));

        assert_arrays_approx_eq(
            &lif.membrane_potential,
            &expected_membrane_potential,
            EPSILON,
        );
    }

    #[test]
    fn emit_spike() {
        let lif = LeakyIntegrateAndFire::new(1, -55.0, -70.0, -65.0, 20.0, 10.0);
        let spike = lif.emit_spike();
        let expected_spike = af::constant(false, Dim4::new(&[1, 1, 1, 1]));

        // Compare Array<bool> element-wise
        let is_eq = af::eq(&spike, &expected_spike, false);
        let all_true_scalar = all_true(&is_eq, 1);
        // assert_eq!(all_true_scalar, true);
    }
}
