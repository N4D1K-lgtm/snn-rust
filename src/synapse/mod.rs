pub mod spike_timing_dependent_plasticity;
pub mod voltage_dependent_plasticity;
pub mod static_synapse;
pub mod connect_synapse;

pub use spike_timing_dependent_plasticity::STDP;
pub use voltage_dependent_plasticity::VoltageDependentPlasticity;
pub use static_synapse::StaticSynapse;
pub use connect_synapse::*;

pub trait Synapse {
    fn update_weight(&mut self, pre_spike: bool, post_spike: bool, dt: f64);
    fn transmit_spike(&self) -> f64;
}