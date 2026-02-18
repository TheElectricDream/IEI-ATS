# IEI-ATS: Inter-Event-Interval Adaptive Time Surface

> An adaptive local time surface algorithm for neuromorphic event camera data, featuring coherence-based noise filtering and temporal surface accumulation in MATLAB.

---

## Overview

The **Inter-Event-Interval Adaptive Time Surface (IEI-ATS)** is a signal processing framework for event-based vision sensors (event cameras / DVS). Unlike conventional frame-based cameras, event cameras asynchronously report per-pixel brightness changes, producing a sparse stream of events `(x, y, t, p)`. IEI-ATS converts this stream into a perceptually coherent, noise-suppressed surface representation suitable for downstream computer vision tasks.

The core contribution is a **three-rule coherence filtering framework** that separates real edge events from noise *before* temporal surface accumulation, combined with a **locally adaptive EMA accumulator** driven by per-pixel inter-event-interval (IEI) statistics. The result is a surface that preserves edge sharpness and polarity contrast while eliminating hot pixels, trailing artifacts, and noise-induced residuals.

---

## Background

### Event Cameras

Event cameras such as the Dynamic Vision Sensor (DVS) and DAVIS sensor output a stream of events rather than frames. Each event encodes:

- `(x, y)` — pixel location
- `t` — timestamp (microsecond resolution)
- `p` — polarity (`+1` for brightness increase, `-1` for decrease)

This asynchronous representation offers high temporal resolution (~1 µs), high dynamic range (~120 dB), and low latency, but requires fundamentally different processing pipelines compared to frame-based methods.

### Time Surfaces

A **time surface** is a 2D map where each pixel stores a decayed representation of when it last fired. The classical exponential time surface takes the form:

$$\mathcal{T}(x, y, t) = \exp\left(-\frac{t - t_{\text{last}}(x,y)}{\tau}\right)$$

where $\tau$ is a fixed time constant controlling the decay rate. This formulation is described in Lagorce et al. (2017), *"HOTS: A Hierarchy of Event-Based Time-Surfaces for Pattern Recognition,"* IEEE TPAMI. IEI-ATS replaces the fixed $\tau$ with a **per-pixel adaptive time constant** derived from local IEI statistics, and replaces the hard reset at each event with an exponential moving average (EMA) update rule.

---

## Algorithm Architecture

The processing pipeline operates on fixed-duration time windows (frames) accumulated from the event stream. Within each frame the following stages are executed in order:

```
Raw Events
    │
    ▼
┌─────────────────────────┐
│   Event Preparation     │  Sort, index, group by pixel
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  Neighborhood Statistics│  IEI mean, std, mean-diff per pixel
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  Coherence Filtering    │  Spatial density · Persistence · IEI regularity
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  Adaptive Time Surface  │  EMA accumulation + divisive normalization
└────────────┬────────────┘
             │
             ▼
    Normalized Output Frame
```

### Stage 1 — Neighborhood Statistics

For each active pixel, the mean and standard deviation of the inter-event interval (IEI) within the current time window are computed. A persistent IEI map is maintained across frames using an **exponential moving average**:

$$\hat{\text{IEI}}_t = (1 - \alpha)\,\hat{\text{IEI}}_{t-1} + \alpha\,\text{IEI}_t$$

where $\alpha$ is the IEI smoothing factor (`iei_alpha`). This retains temporal statistics even in sparse frames.

### Stage 2 — Coherence Filtering

Three independent rules are combined into a scalar coherence score per pixel:

#### Rule 1: Spatial Density (Trace Map)

Events are tested for spatial clustering using a normalized k-d tree radius search over `(x_norm, y_norm, t_norm)` space. The summed distance to spatial neighbours within radius $r_s$ forms a density proxy. Pixels below `trace_threshold` are rejected. The result is log-normalized to compress the heavy-tailed distribution.

#### Rule 2: Temporal Persistence

Active regions are expected to persist across consecutive frames. For each active pixel in the current frame, the nearest active pixel in the previous frame's trace map is found using a 3D KD-tree search over normalized `(row, col, value)` space (`findPersistenceVectorized`). Pixels whose nearest predecessor exceeds `persistence_threshold` are rejected. This eliminates transient noise bursts that do not repeat across windows.

#### Rule 3: IEI Regularity (Similarity Map)

Real edges sweeping across the sensor produce spatially coherent firing rates — adjacent pixels along an edge fire at similar inter-event intervals. Noise events fire at uncorrelated rates relative to their neighbors. The **Coefficient of Variation** (CV = σ/μ) of the local IEI neighborhood is computed via normalized convolution:

$$\text{CV}(x, y) = \frac{\sigma_\text{local}(x, y)}{\mu_\text{local}(x, y)}$$

The regularity score maps CV to `[0, 1]` via:

$$s_\text{reg} = \frac{1}{1 + \text{CV}}$$

An IEI magnitude penalty is then applied to prevent pure-noise regions (where all pixels fire slowly but uniformly) from producing artificially low CV:

$$s_\text{mag} = \exp\left(-\frac{\text{IEI}}{2 \cdot \tilde{\text{IEI}}}\right)$$

where $\tilde{\text{IEI}}$ is the median observed IEI. The combined score is $s = s_\text{reg} \cdot s_\text{mag}$.

The normalized convolution approach follows Knutsson & Westin (1993), *"Normalized and Differential Convolution,"* Proc. IEEE CVPR, pp. 515–523.

The three scores are summed and thresholded at `coherence_threshold` to produce the final binary/soft filter mask.

### Stage 3 — Adaptive Local Time Surface (ALTS)

The ALTS accumulator applies an **asymmetric attack-release envelope** driven by the local IEI statistics:

**Time constant mapping.** The per-pixel IEI is mapped to an effective decay time constant $\tau$ via a sigmoid remap:

$$\tau_\text{active}(x, y) = \sigma_\text{remap}\!\left(\text{IEI}(x,y),\ \tau_\text{min},\ \tau_\text{max}\right)$$

Fast-firing pixels (small IEI) receive small $\tau$ (fast decay); slow-firing pixels receive large $\tau$ (long memory).

**Attack-release envelope.** Active pixels (passing the coherence filter) use $\tau_\text{active}$; inactive pixels with residual surface energy use the fixed release constant $\tau_\text{release} \gg \tau_\text{active}$. This eliminates trailing artifacts from fast-moving edges.

**EMA update (first-order IIR).** The blending coefficient is computed as:

$$\alpha_\text{eff} = 1 - \exp\left(-\frac{\Delta t}{\tau_\text{eff}}\right)$$

This coupling of input gain with decay rate prevents runaway accumulation. The surface update is:

$$S_t = \alpha_\text{eff} \cdot u_t + (1 - \alpha_\text{eff}) \cdot S_{t-1}$$

where $u_t$ is the signed polarity input at the current frame. Active pixels with $\alpha_\text{eff} = 0$ implement a **sample-and-hold**: the previous value is preserved exactly without decay.

**Divisive normalization.** The raw EMA surface is normalized using a variant of the canonical divisive normalization model from Carandini & Heeger (2012), *"Normalization as a canonical neural computation,"* Nature Reviews Neuroscience, 13, 51–62:

$$\hat{S}(x,y) = \frac{S(x,y)}{\sigma + C_\text{smooth}(x,y)^n}$$

where $C_\text{smooth}$ is a spatially pooled event count map, $\sigma$ is the median activity level, and $n = 1$. This compresses the dynamic range without the feedback-loop artifacts that arise when normalized output is fed back into the EMA state.

**Tone mapping.** The final output is rendered using an edge-aware base-detail decomposition. The base layer is compressed via $\tanh$; the detail layer is recombined at a separately controllable boost factor. This prevents flicker artifacts that arise from whole-image normalization.

---

## Repository Structure

```
IEI-ATS/
│
├── main.m                          # Top-level processing script
│
├── +accumulator/
│   ├── localAdaptiveTimeSurface.m  # Core ALTS accumulator (EMA + divisive norm)
│   ├── adaptiveGlobalDecay.m       # Nunes et al. (2023) AGD reference impl.
│   ├── timeSurface.m               # Classical exponential time surface
│   └── speedInvariantTimeSurface.m # Speed-invariant TS (Manderscheid 2019)
│
├── +coherence/
│   ├── computeCoherenceMask.m      # Main coherence pipeline (3-rule combiner)
│   ├── findSpatialNeighbours.m     # KD-tree radius search for spatial density
│   ├── findSimilarities.m          # CV-based IEI regularity scoring
│   ├── findSimilaritiesLOF.m       # LOF-style similarity (experimental)
│   └── findPersistenceVectorized.m # Cross-frame persistence via KNN search
│
├── +stats/
│   ├── computeNeighborhoodStats.m  # Per-pixel IEI mean, std, diff stats
│   ├── spreadEventsSpatially.m     # Spatial splatting for stat computation
│   ├── splatTimestampMap.m         # Gaussian normalized convolution splatting
│   └── printPercentComplete.m      # Progress reporting utility
│
├── +process/
│   ├── sliceToValidRange.m         # Time-window event extraction
│   ├── sigmoidRemap.m              # Sigmoid mapping for tau remapping
│   ├── linearRemap.m               # Linear mapping utility
│   ├── unwrapMap.m                 # Extract filtered events from coherence map
│   └── generateMeshFromFrame.m     # 2D map to point cloud conversion
│
├── +plot/
│   ├── initializeEventVideos.m     # Video writer initialization
│   ├── mapToScatterPlot.m          # 3D scatter visualization
│   ├── mapToSurfPlot.m             # Interpolated surface visualization
│   ├── mapToPCViewer.m             # Point cloud viewer
│   ├── mapToScaledImage.m          # 2D imagesc visualization
│   ├── vectorsToScatterPlot.m      # Scatter from raw vectors
│   └── visualizePersistance.m      # 3D persistence neighbourhood viz
│
└── +voxelization/
    └── discretizeEventsToVoxels.m  # Spatiotemporal voxel occupancy grid
```

---

## Dependencies

- **MATLAB R2021b or later** (uses `imbilatfilt`, `imdilate`, `imgaussfilt`)
- **Image Processing Toolbox** — required for morphological and filter operations
- **Statistics and Machine Learning Toolbox** — required for `createns`, `knnsearch`, `rangesearch`
- Event data in **HDF5 format** with datasets: `/timestamp`, `/x`, `/y`, `/polarity`

To convert AEDAT4 recordings to HDF5, see the companion import script from the [NEXUS](https://github.com/Carleton-SRL/NEXUS/tree/main/import) repository.

---

## Getting Started

### 1. Configure the data path and dataset

In `main.m`, set the path and filename for your HDF5 event recording:

```matlab
hdf5_path = '/path/to/your/datasets/';
file_name  = 'your_recording.hdf5';
```

### 2. Set the processing window

```matlab
t_start_process = 0;    % [seconds]
t_end_process   = 10;   % [seconds]
t_interval      = 0.033; % Accumulation window [seconds]
```

### 3. Tune coherence parameters

The coherence filter parameters are the most dataset-dependent settings:

| Parameter | Description | Typical Range |
|---|---|---|
| `coh_params.r_s` | Spatial density search radius (normalized) | `0.02 – 0.08` |
| `coh_params.trace_threshold` | Min neighbour distance sum to pass density rule | `0.5 – 3.0` |
| `coh_params.persistence_threshold` | Max cross-frame distance to pass persistence rule | `0.0001 – 0.001` |
| `coh_params.coherence_threshold` | Combined score threshold for filter mask | `0.03 – 0.15` |

### 4. Tune time surface parameters

| Parameter | Description | Typical Range |
|---|---|---|
| `alts_params.surface_tau_min` | Minimum decay constant [s] | `0.05 – 0.5` |
| `alts_params.surface_tau_max` | Maximum decay constant [s] | `0.5 – 5.0` |
| `alts_params.surface_tau_release` | Inactive pixel release constant [s] | `1.0 – 5.0` |
| `iei_alpha` | IEI map EMA smoothing factor | `0.5 – 0.95` |

### 5. Run

```matlab
main
```

Output frames are written to `videoWriters{1}` (AVI) and optionally as individual PNGs to the configured `frameOutputFolder`.

---

## Reference Implementations

IEI-ATS includes reference implementations of several published time surface methods for benchmarking:

- **Classical Time Surface** — Lagorce et al. (2017), *"HOTS: A Hierarchy of Event-Based Time-Surfaces for Pattern Recognition,"* IEEE TPAMI, 39(7), 1346–1359.
- **Speed-Invariant Time Surface** — Manderscheid et al. (2019), *"Speed Invariant Time Surface for Learning to Detect Corner Points with Event-Based Cameras,"* Proc. IEEE CVPR. [arXiv:1903.11332](https://arxiv.org/abs/1903.11332)
- **Adaptive Global Decay (AGD)** — Nunes et al. (2023), *"Adaptive Global Decay Process for Event Cameras,"* Proc. IEEE CVPR. [IEEE](https://ieeexplore.ieee.org/document/10205486/)

---

## License

MIT License. Copyright (c) 2026 Alexander Crain. See [LICENSE](LICENSE) for details.

---

## Citation

If you use IEI-ATS in your research, please cite:

```bibtex
@software{crain2026ieats,
  author  = {Crain, Alexander},
  title   = {{IEI-ATS}: Inter-Event-Interval Adaptive Time Surface},
  year    = {2026},
  url     = {https://github.com/[your-github]/IEI-ATS},
  license = {MIT}
}
```
