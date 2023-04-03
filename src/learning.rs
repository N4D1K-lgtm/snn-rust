extern crate arrayfire as af;

use af::{Array, exp, select, gt};


pub trait LearningRule {
    fn update_weights(&mut self, pre_spike_time: &Array<f64>, post_spike_time: &Array<f64>, dt: f64) -> Array<f64>;
}

pub struct STDP {
    pub a_plus: f64,
    pub a_minus: f64,
    pub tau_plus: f64,
    pub tau_minus: f64,
}

impl STDP {
    pub fn new(a_plus: f64, a_minus: f64, tau_plus: f64, tau_minus: f64) -> Self {
        STDP {
            a_plus,
            a_minus,
            tau_plus,
            tau_minus,
        }
    }
}

impl LearningRule for STDP {
    fn update_weights(&mut self, pre_spike_time: &Array<f64>, post_spike_time: &Array<f64>, dt: f64) -> Array<f64> {
        let delta_t = post_spike_time - pre_spike_time;
        let positive_dt = gt(&delta_t, &af::constant(0.0, delta_t.dims()), false);

        let delta_t_dims = delta_t.dims(); 
        let positive_term = af::constant(self.a_plus, delta_t_dims) * exp(&(-delta_t.clone() / self.tau_plus)) * dt;
        let negative_term = af::constant(-self.a_minus, delta_t_dims) * exp(&(-delta_t / self.tau_minus)) * dt;
        
        select(&positive_term, &positive_dt, &negative_term)
    }
}