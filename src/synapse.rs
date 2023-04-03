pub trait Synapse {
    fn update_weight(&mut self, pre_spike: bool, post_spike: bool, dt: f64);
    fn transmit_spike(&self) -> f64;
}

// STDP - Spike-Timing-Dependent Plasticity
pub struct SpikeTimingDependentPlasticity {
    pub weight: f64,
    pub learning_rate: f64,
    pub tau_plus: f64,
    pub tau_minus: f64,
}

impl Synapse for SpikeTimingDependentPlasticity {
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



// StaticSynapse - a synapse with a constant weight
#[derive(Clone)]
pub struct StaticSynapse {
    pub weight: f64,
}

impl Synapse for StaticSynapse {
    fn update_weight(&mut self, _pre_spike: bool, _post_spike: bool, _dt: f64) {
        // The weight does not change in a static synapse
    }

    fn transmit_spike(&self) -> f64 {
        self.weight
    }
}


// VoltageDependentPlasticity
pub struct VoltageDependentPlasticity {
    pub weight: f64,
    pub learning_rate: f64,
    pub voltage_threshold: f64,
}

impl Synapse for VoltageDependentPlasticity {
    fn update_weight(&mut self, pre_spike: bool, post_spike: bool, _dt: f64) {
        if pre_spike && post_spike {
            self.weight +=
                self.learning_rate * (1.0 - (self.weight >= self.voltage_threshold) as i32 as f64);
        }
    }

    fn transmit_spike(&self) -> f64 {
        self.weight
    }
}