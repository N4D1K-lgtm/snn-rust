use super::Synapse;
use crate::learning::{LearningRule, STDP};

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
