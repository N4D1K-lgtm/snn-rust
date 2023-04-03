pub trait Neuron {
    fn new() -> Self;
    fn update_state(&mut self, input_current: f64, dt: f64);
    fn handle_spike(&mut self);
    fn emit_spike(&self) -> bool;
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


/// Hodgkin-Huxley Neuron Model
///
/// Reference: Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of membrane current
/// and its application to conduction and excitation in nerve. The Journal of Physiology, 117(4), 500-544.
///
/// The Hodgkin-Huxley model is a biophysically detailed model that
/// describes the conductance of specific ion channels (sodium and
/// potassium) and their activation and inactivation variables (n, m, and h)
/// to simulate the neuron's behavior more accurately.
///
/// Pros:
/// - High biological realism
/// - Reproduces complex spike patterns
///
/// Cons:
/// - Computationally expensive
/// - Not suitable for large-scale network simulations
///
/// Mathematical equations:
/// dv/dt = (g_na * m^3 * h * (e_na - v) + g_k * n^4 * (e_k - v) + g_l * (e_l - v) + I) / C_m
/// dn/dt = alpha_n * (1 - n) - beta_n * n
/// dm/dt = alpha_m * (1 - m) - beta_m * m
/// dh/dt = alpha_h * (1 - h) - beta_h * h

pub struct HodgkinHuxley {
    pub membrane_potential: f64, // The neuron's membrane potential (voltage) in millivolts (mV)
    pub n: f64, // The activation variable for the delayed rectifier potassium channel (K+)
    pub m: f64, // The activation variable for the fast-activating sodium channel (Na+)
    pub h: f64, // The inactivation variable for the fast-activating sodium channel (Na+)
    pub g_na: f64, // The maximum conductance for the fast-activating sodium channel (Na+) in millisiemens per cm² (mS/cm²)
    pub g_k: f64, // The maximum conductance for the delayed rectifier potassium channel (K+) in millisiemens per cm² (mS/cm²)
    pub g_l: f64, // The leak conductance in millisiemens per cm² (mS/cm²)
    pub e_na: f64, // The sodium Nernst (reversal) potential in millivolts (mV)
    pub e_k: f64, // The potassium Nernst (reversal) potential in millivolts (mV)
    pub e_l: f64, // The leak Nernst (reversal) potential in millivolts (mV)
}

impl HodgkinHuxley {
    pub fn new() -> Self {
        HodgkinHuxley {
            membrane_potential: -65.0,
            n: 0.32,
            m: 0.05,
            h: 0.6,
            g_na: 120.0,
            g_k: 36.0,
            g_l: 0.3,
            e_na: 50.0,
            e_k: -77.0,
            e_l: -54.4,
        }
    }
    fn alpha_n(&self) -> f64 {
        (0.01 * (55.0 + self.membrane_potential))
            / (1.0 - (-0.1 * (55.0 + self.membrane_potential)).exp())
    }

    fn beta_n(&self) -> f64 {
        0.125 * (-65.0 - self.membrane_potential).exp()
    }

    fn alpha_m(&self) -> f64 {
        (0.1 * (40.0 + self.membrane_potential))
            / (1.0 - (-0.1 * (40.0 + self.membrane_potential)).exp())
    }

    fn beta_m(&self) -> f64 {
        4.0 * (-20.0 - self.membrane_potential).exp()
    }

    fn alpha_h(&self) -> f64 {
        0.07 * (-20.0 - self.membrane_potential).exp()
    }

    fn beta_h(&self) -> f64 {
        1.0 / (1.0 + (-0.1 * (30.0 + self.membrane_potential)).exp())
    }
}

impl Neuron for HodgkinHuxley {
    fn new() -> Self {
        Self::new()
    }
    fn update_state(&mut self, input_current: f64, dt: f64) {
        let i_na = self.g_na * self.m.powi(3) * self.h * (self.e_na - self.membrane_potential);
        let i_k = self.g_k * self.n.powi(4) * (self.e_k - self.membrane_potential);
        let i_l = self.g_l * (self.e_l - self.membrane_potential);

        let d_v = (input_current - i_na - i_k - i_l) * dt;

        let dn = (self.alpha_n() * (1.0 - self.n) - self.beta_n() * self.n) * dt;
        let dm = (self.alpha_m() * (1.0 - self.m) - self.beta_m() * self.m) * dt;
        let dh = (self.alpha_h() * (1.0 - self.h) - self.beta_h() * self.h) * dt;

        self.membrane_potential += d_v;
        self.n += dn;
        self.m += dm;
        self.h += dh;
    }

    fn handle_spike(&mut self) {
        self.membrane_potential = -65.0;
        self.n = 0.32;
        self.m = 0.05;
        self.h = 0.6;
    }

    fn emit_spike(&self) -> bool {
        self.membrane_potential >= 30.0
    }
}