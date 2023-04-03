use super::Neuron;

/// Izhikevich Neuron Model
///
/// Reference: Izhikevich, E. M. (2003). Simple model of spiking neurons.
/// IEEE Transactions on Neural Networks, 14(6), 1569-1572.
///
/// The Izhikevich model is a simplified neuron model that can reproduce
/// various spiking patterns observed in biological neurons with a lower
/// computational cost compared to the Hodgkin-Huxley model. The model is
/// governed by two differential equations that describe the dynamics of
/// the membrane potential and a recovery variable.
///
/// Pros:
/// - Reproduces complex spike patterns
/// - Computationally efficient
///
/// Cons:
/// - Limited biological realism
///
/// Mathematical equations:
/// dv/dt = 0.04 * v^2 + 5 * v + 140 - u + I
/// du/dt = a * (b * v - u)
///
/// If v >= 30:
///     v = c
///     u = u + d
///     Emit spike

pub struct Izhikevich {
    pub membrane_potential: f64, // The neuron's membrane potential (voltage) in millivolts (mV)
    pub recovery_variable: f64, // The recovery variable (u) representing the membrane recovery rate
    pub a: f64, // The time scale of the recovery variable (u)
    pub b: f64, // The sensitivity of the recovery variable (u) to the subthreshold fluctuations of the membrane potential
    pub c: f64, // The after-spike reset value of the membrane potential in mV
    pub d: f64, // The after-spike reset value of the recovery variable (u)
}

impl Izhikevich {
    pub fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Izhikevich {
            membrane_potential: -65.0,
            recovery_variable: 0.0,
            a,
            b,
            c,
            d,
        }
    }
}

impl Neuron for Izhikevich {
    fn new() -> Self {
        Self::new(0.02, 0.2, -65.0, 8.0)
    }
    fn update_state(&mut self, input_current: f64, dt: f64) {
        let dv =
            (0.04 * self.membrane_potential.powi(2) +
                5.0 * self.membrane_potential +
                140.0 -
                self.recovery_variable +
                input_current) *
            dt;
        let du = self.a * (self.b * self.membrane_potential - self.recovery_variable) * dt;

        self.membrane_potential += dv;
        self.recovery_variable += du;
    }

    fn handle_spike(&mut self) {
        self.membrane_potential = self.c;
        self.recovery_variable += self.d;
    }

    fn emit_spike(&self) -> bool {
        self.membrane_potential >= 30.0
    }
}