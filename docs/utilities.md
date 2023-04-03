# `utility.rs`

In the `utility.rs` module, you will include utility functions for common tasks, such as generating random numbers, initializing neuron and synapse parameters, and validating input data. You will also implement helper functions for handling common connection patterns (e.g., all-to-all, one-to-one, random) and weight initialization methods (e.g., uniform, Gaussian, constant).

Example:

rust

```rust
pub fn generate_random_numbers(n: usize, min: f64, max: f64) -> Vec<f64> {
    // Generate n random numbers between min and max
}

pub fn validate_connection_pattern(pattern: &str) -> Result<(), String> {
    // Validate the connection pattern string and return an error if it's not supported
}

// Other utility functions and helper functions for initialization and validation
```

These utility functions can be used throughout the library to perform common tasks, reducing code duplication and improving maintainability.