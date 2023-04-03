pub mod leaky_integrate_and_fire;
pub mod hodgkin_huxley;
pub mod izhikevich;

pub use leaky_integrate_and_fire::LeakyIntegrateAndFire;
pub use hodgkin_huxley::HodgkinHuxley;
pub use izhikevich::Izhikevich;
pub trait Neuron {
    fn new() -> Self;
    fn update_state(&mut self, input_current: f64, dt: f64);
    fn handle_spike(&mut self);
    fn emit_spike(&self) -> bool;
}