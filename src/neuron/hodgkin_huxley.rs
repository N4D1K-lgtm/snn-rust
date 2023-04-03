use super::Neuron;

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
