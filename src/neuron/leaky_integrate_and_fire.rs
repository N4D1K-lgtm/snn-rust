use super::Neuron;

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
    pub membrane_potential: f64, // The neuron's membrane potential (voltage) in millivolts (mV)
    pub threshold: f64, // The membrane potential threshold for spike generation in mV
    pub rest_potential: f64, // The resting membrane potential in mV
    pub reset_potential: f64, // The membrane potential reset value after a spike in mV
    pub membrane_resistance: f64, // The neuron's membrane resistance in megaohms (MΩ)
    pub membrane_time_constant: f64, // The membrane time constant in milliseconds (ms), representing the time it takes for the membrane potential to reach ~63% of its final value in response to a step input current
}

impl LeakyIntegrateAndFire {
    pub fn new() -> Self {
        LeakyIntegrateAndFire {
            membrane_potential: -65.0,
            threshold: -50.0,
            rest_potential: -65.0,
            reset_potential: -70.0,
            membrane_resistance: 100.0,
            membrane_time_constant: 10.0,
        }
    }
}

impl Neuron for LeakyIntegrateAndFire {
    fn new() -> Self {
        Self::new()
    }
    fn update_state(&mut self, input_current: f64, dt: f64) {
        let membrane_potential_change =
            ((-(self.membrane_potential - self.rest_potential) +
                self.membrane_resistance * input_current) /
                self.membrane_time_constant) *
            dt;
        self.membrane_potential += membrane_potential_change;
    }

    fn handle_spike(&mut self) {
        self.membrane_potential = self.reset_potential;
    }

    fn emit_spike(&self) -> bool {
        self.membrane_potential >= self.threshold
    }
}