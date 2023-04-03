use crate::neuron::Neuron;

pub struct NeuronPopulation<T: Neuron> {
    pub neurons: Vec<T>,
}

impl<T: Neuron> NeuronPopulation<T> {
    pub fn new(size: usize) -> Self {
        let neurons = (0..size).map(|_| T::new()).collect();
        NeuronPopulation { neurons }
    }

    pub fn update_state(&mut self, input_currents: &[f64], dt: f64) {
        assert_eq!(
            self.neurons.len(),
            input_currents.len(),
            "Number of input currents must match the population size."
        );

        for (neuron, &input_current) in self.neurons.iter_mut().zip(input_currents) {
            neuron.update_state(input_current, dt);
            if neuron.emit_spike() {
                neuron.handle_spike();
            }
        }
    }
}