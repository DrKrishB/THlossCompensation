# THlossCompensation User Manual

**Overview:**
This repository contains MATLAB scripts for analyzing tail neuron hyperactivity patterns and refining motion registration in 2-photon imaging datasets. Core functionality revolves around processing calcium imaging data, tracking tail movement, and applying hierarchical classification logic.

**Requirements:**
MATLAB R2023b or later, Image Processing Toolbox, Signal Processing Toolbox, Suite2p/CaImAn for extraction of nuclear-labeled neurons.

### ⚡ GPU Acceleration

Several scripts in this repository support CUDA-based acceleration using MATLAB’s **Parallel Computing Toolbox**. To take advantage of GPU acceleration:

- Ensure a compatible **NVIDIA GPU** is available on your system.
- Install the latest **CUDA drivers** from NVIDIA and verify that MATLAB detects the GPU using `gpuDevice`.
- Scripts are designed to automatically detect GPU availability and offload computation where possible.

### 📂 File Descriptions (sorted by pipeline order)
| Script                          | Description                                                                                                                       |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `track240fps.m`                 | Extracts tail-tip coordinates from high-frame-rate (240 fps) behavioral videos using adaptive thresholding and centroid tracking. |
| `tracktail2p.m`                 | Aligns extracted tail motion data with 2-photon imaging frames for head-restrained larval zebrafish.                              |
| `suite2pRefineReg.m`           | Refines Suite2p or CaImAn-extracted neuron data, performs anatomical registration via ANTsPy, and maps neural clusters to ~2158 brain regions for quantifying ΔF/F, firing rates, and cell counts. |
| `tailNeuron2x.m`               | Integrates neural activity with tail motion data to assess neuron-behavior coupling.                                               |
| `double_hyperactivity_neurons3D.m` | Identifies hyperactive neurons across 3D space using calcium traces and spatial ROI mapping.                                          |
| `hierarchicalKB.m`             | Applies hierarchical clustering to classify neurons based on extracted activity metrics.                                          |

