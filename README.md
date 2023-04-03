# snn-rust
This Rust library provides tools for developing and simulating spiking neural networks (SNNs). It includes a variety of neuron models and synapse models, as well as tools for creating and managing populations of neurons.

## Installation
To use this library in your Rust project, simply add the following line to your Cargo.toml file:

```makefile
[dependencies]
snn_rust = "0.1.0"
```

## Usage
Here's an example of how to create a population of leaky integrate-and-fire neurons and simulate them:

```rust
// Import necessary modules from your SNN library
use snn_rust::{
    neuron::{LeakyIntegrateAndFire, Neuron},
    network::Network,
    synapse::{SpikeTimingDependentPlasticity, Synapse},
    simulation::Simulation,
};

fn main() {
    // Create neuron populations
    let input_neurons = (0..10)
        .map(|_| LeakyIntegrateAndFire::new(1.0, 1.0, -65.0, -55.0))
        .collect::<Vec<_>>();
    let output_neurons = (0..2)
        .map(|_| LeakyIntegrateAndFire::new(1.0, 1.0, -65.0, -55.0))
        .collect::<Vec<_>>();

    // Create the network and add neuron populations
    let mut network = Network::new();
    let input_pop = network.add_neuron_population(input_neurons);
    let output_pop = network.add_neuron_population(output_neurons);

    // Create synapses with STDP learning rule
    let synapses = (0..input_pop.len() * output_pop.len())
        .map(|_| {
            SpikeTimingDependentPlasticity::new(0.01, 0.01, 0.1, 20.0, 20.0)
        })
        .collect::<Vec<_>>();

    // Connect neuron populations with synapses
    network.connect_populations(
        &input_pop,
        &output_pop,
        synapses,
        crate::snn_lib::utils::all_to_all,
    );

    // Create and configure the simulation
    let mut simulation = Simulation::new(network, 0.1, 1000.0);

    // Run the simulation
    simulation.run();

    // Retrieve and analyze the results
    // E.g., print spike times, neuron states, or synaptic weights
}
```

## Project Structure

```bash
snn_rust/
|-- src/
|   |-- lib.rs
|   |-- neuron.rs
|   |-- synapse.rs
|   |-- network.rs
|   |-- learning.rs
|   |-- simulation.rs
|   |-- utility.rs
|   |-- io.rs
|-- Cargo.toml
```
- `lib.rs`: This is the main entry point for the library. It exports the public API and includes the necessary modules.

- `neuron.rs`: This file contains the Neuron trait definition and implementations of various neuron models (e.g., LeakyIntegrateAndFire, HodgkinHuxley, Izhikevich). Each neuron model is implemented as a separate struct that implements the Neuron trait. The file also includes helper functions for creating neuron populations and initializing their state.

- `synapse.rs`: This file contains the Synapse trait definition and implementations of various synapse models (e.g., StaticSynapse, STDP, VoltageDependentPlasticity). Each synapse model is implemented as a separate struct that implements the Synapse trait. The file also includes helper functions for creating and managing synaptic connections between neuron populations.

- `network.rs`: This file contains the Network struct, which represents the SNN and stores neuron populations, synaptic connections, and additional network properties. It includes methods for adding neuron populations, creating synaptic connections between populations, and updating the network state at each time step.

- `learning.rs`: This file contains the LearningRule trait definition and implementations of various learning rules (e.g., STDP, RSTDP, HomeostaticPlasticity). Each learning rule is implemented as a separate struct that implements the LearningRule trait. This file is imported by the synapse.rs module to allow users to specify the learning rule when creating synaptic connections.

- `simulation.rs`: This file contains the Simulation struct, which manages the simulation of the SNN, including the network state, simulation parameters, and time step. It includes methods for initializing the simulation, updating the network state at each time step, and recording the simulation results.

- `utility.rs`: This file contains utility functions for common tasks, such as generating random numbers, initializing neuron and synapse parameters, and validating input data. It also includes helper functions for handling common connection patterns and weight initialization methods.

- `io.rs`: This file contains functions for loading and saving network configurations, neuron and synapse parameters, and simulation results. It supports various file formats, such as JSON, CSV, or binary, for easy interoperability with other tools and libraries.

## Module Dependencies:

neuron.rs, synapse.rs, and learning.rs are imported by network.rs to create and manage the SNN.
network.rs is imported by simulation.rs to simulate the SNN.
utility.rs and io.rs are imported by other modules as needed for various utility functions and file I/O operations.

## Crate Dependencies:

- rand = "0.8.5"
- rand_distr = "0.4.3"
- ndarray = "0.15.6"
