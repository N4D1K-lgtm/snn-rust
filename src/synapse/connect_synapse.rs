// Import the necessary modules
use rand_distr::{Distribution, Normal, Uniform};
use rand::Rng;

use super::Synapse;
use crate::neuron::Neuron;


#[derive(Clone)]
pub struct Connection<S: Synapse + Clone> {
    pub synapse_index: usize,
    pub weight: f64,
    pub synapse: S,
}

pub enum ConnectionPattern {
    AllToAll,
    OneToOne,
    Random(f64),
}

pub enum WeightInitialization {
    Constant(f64),
    Uniform(f64, f64),
    Gaussian(f64, f64),
}

pub fn connect_synapse<S: Synapse + Clone>(
    pre_population: &[impl Neuron],
    post_population: &mut [impl Neuron],
    connection_matrix: &[Vec<Connection<S>>],
) {
    for (pre_index, pre_neuron) in pre_population.iter().enumerate() {
        if pre_neuron.emit_spike() {
            for connection in &connection_matrix[pre_index] {
                let input_current = connection.synapse.transmit_spike() * connection.weight;
                post_population[connection.synapse_index].update_state(input_current, 0.0);
            }
        }
    }
}


pub fn create_connections<S: Synapse + Clone>(
    pre_population_size: usize,
    post_population_size: usize,
    synapse_model: S,
    pattern: ConnectionPattern,
    weight_initialization: WeightInitialization,
) -> Vec<Vec<Connection<S>>> {
    let mut rng = rand::thread_rng();

    let mut connection_matrix: Vec<Vec<Connection<S>>> = vec![Vec::new(); pre_population_size];


    match pattern {
        ConnectionPattern::AllToAll => {
            for i in 0..pre_population_size {
                for j in 0..post_population_size {
                    let weight = generate_weight(&mut rng, &weight_initialization);
                    let connection = Connection {
                        synapse_index: j,
                        weight,
                        synapse: synapse_model.clone(),
                    };
                    connection_matrix[i].push(connection);
                }
            }
        }
        ConnectionPattern::OneToOne => {
            for i in 0..pre_population_size.min(post_population_size) {
                let weight = generate_weight(&mut rng, &weight_initialization);
                let connection = Connection {
                    synapse_index: i,
                    weight,
                    synapse: synapse_model.clone(),
                };
                connection_matrix[i].push(connection);
            }
        }
        ConnectionPattern::Random(probability) => {
            for i in 0..pre_population_size {
                for j in 0..post_population_size {
                    if rng.gen::<f64>() < probability {
                        let weight = generate_weight(&mut rng, &weight_initialization);
                        let connection = Connection {
                            synapse_index: j,
                            weight,
                            synapse: synapse_model.clone(),
                        };
                        connection_matrix[i].push(connection);
                    }
                }
            }
        }
    }

    connection_matrix
}

fn generate_weight(rng: &mut impl Rng, weight_initialization: &WeightInitialization) -> f64 {
    match weight_initialization {
        WeightInitialization::Constant(value) => *value,
        WeightInitialization::Uniform(min, max) => {
            let distribution = Uniform::new(*min, *max);
            distribution.sample(rng)
        }
        WeightInitialization::Gaussian(mean, std_dev) => {
            let distribution = Normal::new(*mean, *std_dev)
                .expect("Failed to create Gaussian distribution");
            distribution.sample(rng)
        }
    }
}

// // Create a tests module
// #[cfg(test)]
// mod tests {
//     use super::*;

//     // Test the behavior of StaticSynapse
//     #[test]
//     fn test_static_synapse() {
//         let synapse = StaticSynapse { weight: 0.5 };
//         assert_eq!(synapse.transmit_spike(), 0.5);
//     }

//     // Test the behavior of STDP
//     #[test]
//     fn test_stdp() {
//         let mut synapse = STDP {
//             weight: 0.0,
//             learning_rate: 0.01,
//             tau_plus: 20.0,
//             tau_minus: 20.0,
//         };

//         synapse.update_weight(true, true, 5.0);
//         assert!((synapse.weight - 0.000247009).abs() < 1e-8);

//         synapse.update_weight(true, false, 5.0);
//         assert!((synapse.weight - 0.000494018).abs() < 1e-8);
//     }

//     // Test the behavior of VoltageDependentPlasticity
//     #[test]
//     fn test_voltage_dependent_plasticity() {
//         let mut synapse = VoltageDependentPlasticity {
//             weight: 0.0,
//             learning_rate: 0.01,
//             voltage_threshold: 0.5,
//         };

//         synapse.update_weight(true, true, 0.0);
//         assert_eq!(synapse.weight, 0.01);
//     }

//     // Test connecting neuron populations using StaticSynapse
//     #[test]
//     fn test_connect_synapse_static() {
//         let pre_population: Vec<LeakyIntegrateAndFire> = vec![LeakyIntegrateAndFire::new(); 2];
//         let mut post_population: Vec<LeakyIntegrateAndFire> = vec![LeakyIntegrateAndFire::new(); 2];

//         let synapse_model = StaticSynapse::new(0.5);
//         let pattern = ConnectionPattern::AllToAll;
//         let weight_initialization = WeightInitialization::Constant(0.5);

//         let connection_matrix = create_connections(
//             pre_population.len(),
//             post_population.len(),
//             synapse_model,
//             pattern,
//             weight_initialization,
//         );

//         connect_synapse(&pre_population, &mut post_population, &connection_matrix);

//         // Check the updated membrane potentials of post-synaptic neurons
//         assert_eq!(post_population[0].membrane_potential, 0.0);
//         assert_eq!(post_population[1].membrane_potential, 0.0);
//     }
// }