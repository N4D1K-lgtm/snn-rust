# `neuron.rs`

In the `neuron.rs` module, you will define the `Neuron` trait and implement various spiking neuron models as structs. The `Neuron` trait should outline the common interface for all neuron models, including methods for updating the neuron state, handling incoming spikes, and emitting spikes.

Example:

rust

```rust
pub trait Neuron {
    fn update_state(&mut self, dt: f64);
    fn handle_incoming_spike(&mut self, weight: f64);
    fn emit_spike(&self) -> bool;
}

pub struct LeakyIntegrateAndFire {
    // Parameters and state variables specific to this model
}

impl LeakyIntegrateAndFire {
    // Additional methods specific to this model, if any
}

impl Neuron for LeakyIntegrateAndFire {
    // Implement the Neuron trait methods for LeakyIntegrateAndFire
}

// Implement similar structs and impl blocks for HodgkinHuxley and Izhikevich models.
```

The `neuron.rs` module should also include functions to create neuron populations, initialize their state, and configure their parameters. These functions can be implemented as methods or associated functions of the structs or as standalone functions in the module.