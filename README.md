# snn-rust
This Rust library provides tools for developing and simulating spiking neural networks (SNNs). It includes a variety of neuron models and synapse models, as well as tools for creating and managing populations of neurons.

## Installation
To use this library in your Rust project, simply add the following line to your Cargo.toml file:

```makefile
snn = "0.1.0"
```

## Usage
Here's an example of how to create a population of leaky integrate-and-fire neurons and simulate them:

```rust
use snn::neuron::{LeakyIntegrateAndFire, Neuron};
use snn::network::NeuronPopulation;
use rand::prelude::*;

fn main() {
    let mut rng = rand::thread_rng();
    let mut lif_population = NeuronPopulation::<LeakyIntegrateAndFire>::new(10);

    let input_currents = vec![0.0; 10];

    for i in 0..100 {
        let dt = 0.1;
        let input_currents = input_currents.as_slice();
        lif_population.update_state(input_currents, dt);

        for neuron in &lif_population.neurons {
            if neuron.emit_spike() {
                println!("Neuron {} fired!", i);
            }
        }
    }
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
