numpress
========

[![Build Status](https://api.travis-ci.org/Alexhuszagh/numpress.svg?branch=master)](https://travis-ci.org/Alexhuszagh/numpress)
[![Latest Version](https://img.shields.io/crates/v/numpress.svg)](https://crates.io/crates/numpress)

A pure rust implementation of [ms-numpress](https://github.com/ms-numpress/ms-numpress), a fast, minimally lossy compression algorithm for mass spectrometry data.

# Getting Started

```rust
// Compress floats to bytes
let floats: Vec<f64> = vec![100., 101., 102., 103.];
let compressed: Vec<u8> = numpress_compress(&decoded, DEFAULT_SCALING)?;

// Decompress floats from bytes.
let decompressed: Vec<f64> = numpress_decompress(&compressed)?;
```

# Documentation

Numpress's documentation can be found on [docs.rs](https://docs.rs/numpress).

# Dependency

Numpress is available on crates.io. Use the following in Cargo.toml:

```yaml
[dependencies]
numpress = "1.0"
```

# License

Like the original [ms-numpress](https://github.com/ms-numpress/ms-numpress) implementation, this code is open source. It is dual licensed under the Apache 2.0 license as well as the 3-clause BSD license. See the LICENCE-BSD and the LICENCE-APACHE files for the licenses.

# Contributing

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in numpress by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any additional terms or conditions.
