extern crate arrayfire as af;

use crate::learning::LearningRule;
use af::{Array, Dim4};
pub trait Synapse {
    fn update_weight(&mut self, pre_spike: &Array<bool>, post_spike: &Array<bool>, dt: f64);
    fn transmit_spike(&self) -> Array<f64>;
}

pub struct SpikeTimingDependentPlasticity {
    pub weight: Array<f64>,
    learning_rule: Box<dyn LearningRule>,
    pub pre_spike_time: Array<f64>,
    pub post_spike_time: Array<f64>,
}

impl SpikeTimingDependentPlasticity {
    pub fn new(
        n_synapses: usize,
        initial_weight: f64,
        learning_rule: Box<dyn LearningRule>,
    ) -> Self {
        SpikeTimingDependentPlasticity {
            weight: af::constant(
                initial_weight,
                Dim4::new(&[n_synapses.try_into().unwrap(), 1, 1, 1]),
            ),
            learning_rule,
            pre_spike_time: af::constant(
                f64::NEG_INFINITY,
                Dim4::new(&[n_synapses.try_into().unwrap(), 1, 1, 1]),
            ),
            post_spike_time: af::constant(
                f64::NEG_INFINITY,
                Dim4::new(&[n_synapses.try_into().unwrap(), 1, 1, 1]),
            ),
        }
    }
}

impl Synapse for SpikeTimingDependentPlasticity {
    fn update_weight(&mut self, pre_spike: &Array<bool>, post_spike: &Array<bool>, dt: f64) {
        let spike_occurred = af::or(&pre_spike, &post_spike, false);
        let weight_change =
            self.learning_rule
                .update_weights(&self.pre_spike_time, &self.post_spike_time, dt);
        let updated_weight = &self.weight + weight_change;
        self.weight = af::select(&updated_weight, &spike_occurred, &self.weight);
    }

    fn transmit_spike(&self) -> Array<f64> {
        self.weight.clone()
    }
}

// // StaticSynapse - a synapse with a constant weight
// #[derive(Clone)]
// pub struct StaticSynapse {
//     pub weight: f64,
// }

// impl Synapse for StaticSynapse {
//     fn update_weight(&mut self, _pre_spike: bool, _post_spike: bool, _dt: f64) {
//         // The weight does not change in a static synapse
//     }

//     fn transmit_spike(&self) -> f64 {
//         self.weight
//     }
// }

// // VoltageDependentPlasticity
// pub struct VoltageDependentPlasticity {
//     pub weight: f64,
//     pub learning_rate: f64,
//     pub voltage_threshold: f64,
// }

// impl Synapse for VoltageDependentPlasticity {
//     fn update_weight(&mut self, pre_spike: bool, post_spike: bool, _dt: f64) {
//         if pre_spike && post_spike {
//             self.weight +=
//                 self.learning_rate * (1.0 - (self.weight >= self.voltage_threshold) as i32 as f64);
//         }
//     }

//     fn transmit_spike(&self) -> f64 {
//         self.weight
//     }
// }
