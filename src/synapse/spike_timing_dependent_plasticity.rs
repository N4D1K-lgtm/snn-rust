use super::Synapse;

// STDP - Spike-Timing-Dependent Plasticity
pub struct STDP {
    pub weight: f64,
    pub learning_rate: f64,
    pub tau_plus: f64,
    pub tau_minus: f64,
}

impl Synapse for STDP {
    fn update_weight(&mut self, pre_spike: bool, post_spike: bool, dt: f64) {
        if pre_spike {
            self.weight += self.learning_rate
                * (post_spike as i32 * 2 - 1) as f64
                * f64::exp(-dt / self.tau_plus);
        } else if post_spike {
            self.weight += self.learning_rate
                * (pre_spike as i32 * 2 - 1) as f64
                * f64::exp(-dt / self.tau_minus);
        }
    }

    fn transmit_spike(&self) -> f64 {
        self.weight
    }
}
