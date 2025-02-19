# MonoAlg3D ![build](https://github.com/rsachetto/MonoAlg3D_C/actions/workflows/build.yml/badge.svg)

The MonoAlg3D is a program for solving the 3D monodomain equation by applying the Finite Volume Method.

# Pre-Requisites

  - Linux-based operating system
  - Intel compiler [SYCL version]
  - NVIDIA driver [CUDA version]
  - CUDA compiler [CUDA version]

## Install Intel oneAPI for Fedora 40

- Go to this link and download the Intel oneAPI Base Toolkit for Linux Offline Installer.
- https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?packages=oneapi-toolkit&oneapi-toolkit-os=linux&oneapi-lin=offline

- Create a folder called _"OneAPI-Install"_ in your _/home_ directory and move the file you downloaded to here

```sh
$ cd; mkdir OneAPI-Install
$ mv Downloads/intel-oneapi-base-toolkit-2025.0.1.46_offline.sh OneAPI-Install
```

- Execute the setup using the GUI and follow all the instructions

```sh
$ cd OneAPI-Install
$ sudo sh ./intel-oneapi-base-toolkit-2025.0.1.46_offline.sh
```

## Install the Codeplay NVIDIA compatibility plugin

- This guide contains information on using DPC++ to run SYCL™ applications on NVIDIA® GPUs via the DPC++ CUDA® plugin.

- Go to this website: https://developer.codeplay.com/products/oneapi/nvidia/download/

- Download a version according to your GPU device and CUDA version

- Execute the shell script and install everything

```sh
$ sudo sh oneapi-for-nvidia-gpus-2025.0.0-linux.sh
``` 

- Load the enviroment variables

```sh
$ . /opt/intel/oneapi/setvars.sh --include-intel-llvm
```

- Check if SYCL can identify your GPU now with the command:

```sh
$ sycl-ls
[opencl:cpu][opencl:0] Intel(R) OpenCL, Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz OpenCL 3.0 (Build 0) [2024.18.12.0.05_160000]
[cuda:gpu][cuda:0] NVIDIA CUDA BACKEND, NVIDIA GeForce GTX 1060 6GB 6.1 [CUDA 12.7]
```

## Run examples

- Start the enviroment by running

```sh
$ . /opt/intel/oneapi/setvars.sh --include-intel-llvm
```

### Compile

#### CUDA version

```sh
$ ./build.sh -f
```

#### OpenMP version

```sh
$ ./build.sh -f only_cpu
```

#### SYCL version

```sh
$ ./build.sh -f sycl
```
The binary files will be saved in the ```bin``` folder.

### Setting the enviroment

Ubuntu: Refer to [the ubuntu guide](guide-monoalg3d-ubuntu.md)

Fedora: Refer to [the fedora guide](guide-monoalg3d-fedora.md)

Windows: Refer to [the windows guide](guide-monoalg3d-windows.md)

# Running examples
```sh
$ bin/MonoAlg3D -c example_configs/cuboid_ohara.ini 
```

The output will be saved in the VTK format. In order to see the results you can use Paraview (https://www.paraview.org/) or the compiled visualization tool in ```bin/MonoAlg3D_visualizer```. You can also set the output to plain text, by changing the section ```save_result``` in example_configs/cuboid_ohara.ini to:

```ìni
[save_result]
print_rate=250
output_dir=./outputs/tmp_cube
main_function=save_as_text_or_binary
binary=false
```

In the plain text format we have:

- Each line represents a Volume
- Each volume is represented by its center point (X, Y, and Z), the value of half of its side length on x, y and z and the calculated V

Example file:

```
850,850,950,50,50,50, -85
850,950,950,50,50,50, -85
850,950,850,50,50,50, -85
```

This file represents 3 volumes with 100 micrometer of side. The first volume is centered at  at 850,850,950 and the calculated V is -85 mV.

# Contributors:

@rsachetto Rafael Sachetto Oliveira

@bergolho Lucas Arantes Berg

@Rodrigo-Weber-dos-Santos Rodrigo Weber dos Santos

Among others.

# How to cite:

Oliveira RS, Rocha BM, Burgarelli D, Meira Jr W, Constantinides C, dos Santos RW. Performance evaluation of GPU parallelization, space‐time adaptive algorithms, and their combination for simulating cardiac electrophysiology. Int J Numer Meth Biomed Engng. 2018;34:e2913. https://doi.org/10.1002/cnm.2913

# Credits
[Heart icons created by phatplus - Flaticon](https://www.flaticon.com/free-icons/heart)
