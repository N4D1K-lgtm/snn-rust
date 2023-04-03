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
    fn update_synaptic_weight(&mut self, pre_spike_time: f64, post_spike_time: f64) -> f64 {
        let delta_t = post_spike_time - pre_spike_time;

        if delta_t > 0.0 {
            self.a_plus * (-delta_t / self.tau_plus).exp()
        } else {
            -self.a_minus * (-delta_t / self.tau_minus).exp()
        }
    }
}