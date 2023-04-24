# Parallel_GSP

This repository contains code for implementing the sampling and reconstruction techniques presented in [Parallel Graph Signal Processing: Sampling and Reconstruction](https://ieeexplore.ieee.org/abstract/document/10081083).

## Reconstruction

`bmrp.m` function implements blueness minimization reconstruction with partitions. The inputs are:

* `x`: a noisy signal vector of length N with entries representing the nodes' signal value.
* `sampling_nodes`: a set of sampling nodes.
* `L`: the Laplacian matrix.
* `epsilon`: a perturbation parameter.
* `q`: the Laplacian power.
* `partitions`: a vector of length N with entries representing the nodes' partition.

The output is:

* `x_rec`: the reconstructed signal.

## PVAC

`pvac.m` function implements partitioned void and cluster sampling. It requires the `error_diffusion.mat` and `my_vca_v00.m` files. The inputs are:

* `W`: the weight matrix.
* `n_sample`s: the number of samples.
* `partition`s: a vector of length N with entries representing the nodes' partition.

The output is:
* `sampling_nodes`: a set of sampling nodes.

## DD

`dot_diffusion.m` function implements dot diffusion sampling. It requires the `dot_diffusion_partition.m` file. The inputs are:

* `W`: the weight matrix.
* `n_samples`: the number of samples.
* `partitions`: a vector of length N with entries representing the nodes' partition.

The output is:

* `sampling_nodes`: a set of sampling nodes.

##Citation:

If you use this code in your research, please cite the following paper:

Dapena, Daniela, Daniel L. Lau, and Gonzalo R. Arce. "Parallel Graph Signal Processing: Sampling and Reconstruction." IEEE Transactions on Signal and Information Processing over Networks (2023).
