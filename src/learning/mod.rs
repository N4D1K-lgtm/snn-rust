pub mod 
pub use STDP::STDP;
pub trait LearningRule {
    fn update_synaptic_weight(&mut self, pre_spike_time: f64, post_spike_time: f64) -> f64;
}