# `lib.rs`

In this file, you will declare and export the public API of your library. You will import the necessary modules and re-export their public types, traits, and functions. This will serve as the entry point for users of your library to interact with the components.

Example:

rust

```rust
pub mod neuron;
pub mod synapse;
pub mod network;
pub mod learning;
pub mod simulation;
pub mod utility;
pub mod io;

// Re-export public types, traits, and functions.
pub use neuron::{Neuron, LeakyIntegrateAndFire, HodgkinHuxley, Izhikevich};
pub use synapse::{Synapse, StaticSynapse, STDP, VoltageDependentPlasticity};
pub use network::Network;
pub use learning::{LearningRule, STDP, RSTDP, HomeostaticPlasticity};
pub use simulation::Simulation;
pub use utility::{UtilityFunction1, UtilityFunction2, UtilityFunction3};
pub use io::{load_network, save_network};
```