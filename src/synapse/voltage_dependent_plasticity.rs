use super::Synapse;

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