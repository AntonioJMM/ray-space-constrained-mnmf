# Ray-Space Constrained Multichannel Nonnegative Matrix Factorization for Audio Source Separation

This repository contains the implementation of the **Ray-Space Constrained Multichannel Nonnegative Matrix Factorization (MNMF)** method for **audio source separation**. This approach extends conventional MNMF by integrating a **Ray-Space model**, enabling a more precise representation of spatial characteristics and improving source separation performance.

## Features

- **Ray-Space Dictionary Construction**: Models the propagation of signals using Green’s functions and applies a Ray-Space transform.
- **Frequency-Dependent Propagation Modeling**: Addresses the limitations of traditional far-field models by incorporating frequency dependency.
- **Regularized Source Activation**: Introduces constraints to prevent multiple sources from being assigned to the same grid position simultaneously.

## Installation

Clone this repository and set up the required dependencies:

```
git clone https://github.com/your_username/ray-space-constrained-mnmf.git
cd ray-space-constrained-mnmf
```

Ensure that MATLAB is installed on your system and add the repository path to MATLAB using the following span

## Usage

Run the `main.m` script to perform source separation on an example dataset

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature-name`).
3. Commit your changes (`git commit -m "Add feature"`).
4. Push to your branch (`git push origin feature-name`).
5. Open a pull request.

## Reference

If you use this implementation in your research, please cite:

> A. J. Muñoz-Montoro, M. Olivieri, M. Pezzoli, J. Carabias-Orti, F. Antonacci and A. Sarti, **"Ray-Space Constrained Multichannel Nonnegative Matrix Factorization for Audio Source Separation,"**  *2024 32nd European Signal Processing Conference (EUSIPCO)* , Lyon, France, 2024, pp. 396-400, doi: [10.23919/EUSIPCO63174.2024.10715403](https://doi.org/10.23919/EUSIPCO63174.2024.10715403).

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
