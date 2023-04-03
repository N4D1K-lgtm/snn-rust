use rand::Rng;
use rand_distr::{Distribution, Normal, Uniform};

/// Trait for defining a neuron's behavior
pub trait Neuron {
    /// Updates the state of the neuron given input current and delta time
    fn update_state(&mut self, input_current: f64, dt: f64);

    /// Determines if the neuron emits a spike
    fn emit_spike(&self) -> bool;

    /// Handles the emission of a spike by the neuron
    fn handle_spike(&mut self);
}

/// Trait for defining a synapse's behavior
pub trait Synapse {
    /// Transmits a spike and returns the current induced in the post-synaptic neuron
    fn transmit_spike(&mut self) -> f64;
}

/// A network of neurons and synapses
pub struct Network<T: Neuron, S: Synapse> {
    populations: Vec<NeuronPopulation<T>>,
    connections: Vec<Vec<Connection<S>>>,
}

impl<T: Neuron, S: Synapse> Network<T, S> {
    /// Creates a new empty network
    pub fn new() -> Self {
        Network {
            populations: Vec::new(),
            connections: Vec::new(),
        }
    }

    /// Adds a population of neurons to the network and returns the index of the population
    pub fn add_neuron_population(&mut self, population: NeuronPopulation<T>) -> usize {
        let index = self.populations.len();
        self.populations.push(population);
        index
    }

    /// Connects two populations of neurons with synapses according to the specified pattern and weight initialization
    pub fn connect_populations(
        &mut self,
        pre_population_index: usize,
        post_population_index: usize,
        synapse_model: S,
        connection_pattern: ConnectionPattern,
        weight_initialization: WeightInitialization,
    ) {
        let pre_population = &self.populations[pre_population_index];
        let post_population = &mut self.populations[post_population_index];

        let connection_matrix = create_connections(
            pre_population.neurons.len(),
            post_population.neurons.len(),
            synapse_model,
            connection_pattern,
            weight_initialization,
        );

        self.connections.push(connection_matrix);
    }

    /// Updates the state of all neurons in the network given an array of input currents and delta time
    pub fn update_state(&mut self, input_currents: &[f64], dt: f64) {
        for (i, population) in self.populations.iter_mut().enumerate() {
            population.update_state(input_currents[i], dt);
        }

        for (i, connection_matrix) in self.connections.iter().enumerate() {
            let pre_population_index = i;
            let post_population_index = i;

            let pre_population = &self.populations[pre_population_index];
            let post_population = &mut self.populations[post_population_index];

            connect_synapse(pre_population, post_population, connection_matrix);
        }
    }
}

/// A container for a population of neurons of type `T`.
pub struct NeuronPopulation<T: Neuron> {
    /// A vector containing the neurons in the population.
    neurons: Vec<T>,
}

impl<T: Neuron> NeuronPopulation<T> {
    /// Creates a new `NeuronPopulation` with the specified number of neurons.
    ///
    /// # Arguments
    ///
    /// * `size` - The number of neurons in the population.
    ///
    /// # Returns
    ///
    /// A new `NeuronPopulation` with the specified number of neurons.
    pub fn new(size: usize) -> Self {
        let neurons = (0..size).map(|_| T::new()).collect();
        NeuronPopulation { neurons }
    }

    /// Updates the state of all neurons in the population based on the input current and time step.
    ///
    /// # Arguments
    ///
    /// * `input_current` - The input current to apply to all neurons in the population.
    /// * `dt` - The time step to use for updating the neurons' states.
    pub fn update_state(&mut self, input_current: f64, dt: f64) {
        for neuron in &mut self.neurons {
            neuron.update_state(input_current, dt);
            if neuron.emit_spike() {
                neuron.handle_spike();
            }
        }
    }
}

/// A struct representing a connection between two neurons.
#[derive(Clone)]
pub struct Connection<S: Synapse> {
    /// The index of the post-synaptic neuron.
    synapse_index: usize,
    /// The weight of the connection.
    weight: f64,
    /// The synapse used to transmit spikes across the connection.
    synapse: S,
}

/// An enumeration of possible connection patterns between neuron populations.
pub enum ConnectionPattern {
    /// All-to-all connection pattern.
    AllToAll,
    /// One-to-one connection pattern.
    OneToOne,
    /// Random connection pattern with the specified probability of connection.
    Random(f64),
}

/// An enumeration of possible weight initialization methods.
pub enum WeightInitialization {
    /// All connections have a constant weight.
    Constant(f64),
    /// Connections have weights drawn uniformly from a specified range.
    Uniform(f64, f64),
    /// Connections have weights drawn from a Gaussian distribution with the specified mean and standard deviation.
    Gaussian(f64, f64),
}

/// Connects a pre-synaptic neuron population to a post-synaptic neuron population using the specified connection matrix.
///
/// # Arguments
///
/// * `pre_population` - The pre-synaptic neuron population.
/// * `post_population` - The post-synaptic neuron population.
/// * `connection_matrix` - The connection matrix specifying the connections between neurons.
pub fn connect_synapse<T: Neuron, S: Synapse + Clone>(
    pre_population: &NeuronPopulation<T>,
    post_population: &mut NeuronPopulation<T>,
    connection_matrix: &[Connection<S>],
) {
    for (pre_index, pre_neuron) in pre_population.neurons.iter().enumerate() {
        if pre_neuron.emit_spike() {
            for connection in &connection_matrix[pre_index..pre_index + 1] {
                let input_current = connection.synapse.transmit_spike() * connection.weight;
                post_population
                    .neurons
                    .get_mut(connection.synapse_index)
                    .unwrap()
                    .update_state(input_current, 0.0);
            }
        }
    }
}

/// Creates a matrix of synapse connections between two neuron populations, according to a pattern
/// and with specified weight initialization.
///
/// # Arguments
///
/// * `pre_population_size` - The size of the pre-synaptic neuron population.
/// * `post_population_size` - The size of the post-synaptic neuron population.
/// * `synapse_model` - The synapse model used to create each connection.
/// * `pattern` - The connection pattern between the pre and post populations.
/// * `weight_initialization` - The method of weight initialization.
///
/// # Returns
///
/// A vector of `Connection` structures representing the synapse connections between the two
/// neuron populations.
pub fn create_connections<S: Synapse + Clone>(
    pre_population_size: usize,
    post_population_size: usize,
    synapse_model: S,
    pattern: ConnectionPattern,
    weight_initialization: WeightInitialization,
) -> Vec<Connection<S>> {
    // Generate a random number generator instance
    let mut rng = rand::thread_rng();

    // Create a vector to hold the connections
    let mut connection_matrix: Vec<Connection<S>> =
        Vec::with_capacity(pre_population_size * post_population_size);

    // Generate the connections according to the given pattern
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
                    connection_matrix.push(connection);
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
                connection_matrix.push(connection);
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
                        connection_matrix.push(connection);
                    }
                }
            }
        }
    }

    connection_matrix
}

/// Generate a weight for a synapse connection according to a specified initialization method.
///
/// # Arguments
///
/// * `rng` - A mutable reference to a random number generator instance.
/// * `weight_initialization` - The method of weight initialization.
///
/// # Returns
///
/// A `f64` representing the weight of the synapse connection.
fn generate_weight(rng: &mut impl Rng, weight_initialization: &WeightInitialization) -> f64 {
    match weight_initialization {
        WeightInitialization::Constant(value) => *value,
        WeightInitialization::Uniform(min, max) => {
            let distribution = Uniform::new(*min, *max);
            distribution.sample(rng)
        }
        WeightInitialization::Gaussian(mean, std_dev) => {
            let distribution =
                Normal::new(*mean, *std_dev).expect("Failed to create Gaussian distribution");
            distribution.sample(rng)
        }
    }
}
