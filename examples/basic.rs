use snn_rust::neuron::{LeakyIntegrateAndFire, Neuron};
use snn_rust::network::{NeuronPopulation};

fn main() {
    println!("hello from rust");

    let mut lif_population = NeuronPopulation::<LeakyIntegrateAndFire>::new(10);

    // Increase the input current to 15.0 for each neuron
    let input_currents = vec![15.0; lif_population.neurons.len()];
    let dt = 0.1; // Time step in milliseconds

    // Increase the number of iterations to 1000 for a longer simulation
    for t in 0..1000 {
        lif_population.update_state(&input_currents, dt);
        for (i, neuron) in lif_population.neurons.iter_mut().enumerate() {
            if neuron.emit_spike() {
                // Print the neuron index and simulation time when the neuron spikes
                println!("Neuron {} spiked at time: {} ms", i, t as f64 * dt);
                neuron.handle_spike();
            }
        }
    }
}