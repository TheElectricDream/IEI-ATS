# IEI-ATS: Inter-Event-Interval Adaptive Time Surface

> An adaptive local time surface algorithm for neuromorphic event camera data, featuring coherence-based noise filtering and temporal surface accumulation in MATLAB.

<div align="center">
    <img src="assets/iei-ats-demo.gif" alt="Demo"/>
    <p><strong>Figure 1.</strong> <em>IEI-ATS output surface on the EVOS dataset. Rotating target spacecraft under nominal lighting conditions.</em></p>
</div>

---

> [!IMPORTANT]
> This repository represents the _current best implementation_ of the IEI-ATS algorithm. This means that the algorithm is an improved version of the one described in the related publication (Event-Based Spacecraft Representation Using Inter-Event-Interval Adaptive Time Surfaces). However, the core principle of the approach remains the same.

## Overview

The **Inter-Event-Interval Adaptive Time Surface (IEI-ATS)** is a signal processing framework for event-based vision sensors (event cameras / DVS). Unlike conventional frame-based cameras, event cameras asynchronously report per-pixel brightness changes, producing a sparse stream of events `(x, y, t, p)`. IEI-ATS converts this stream into a perceptually coherent, noise-suppressed surface representation suitable for downstream computer vision tasks.

The core contribution is a **three-rule coherence filtering framework** that separates real edge events from noise *before* temporal surface accumulation, combined with a **locally adaptive EMA accumulator** driven by per-pixel inter-event-interval (IEI) statistics. The result is a surface that preserves edge sharpness and polarity contrast while eliminating hot pixels, trailing artifacts, and noise-induced residuals.

This work also contributes the [EVent-based Observation of Spacecraft (EVOS)](https://www.frdr-dfdr.ca/repo/dataset/1c83d3dc-3a8c-408b-a6c4-9251c0ec1d05) dataset. This dataset is designed to support research in event-based vision for autonomous on-orbit inspection and space debris removal. The dataset was collected at Carleton University using the Spacecraft Proximity Operations Testbed (SPOT). The data was captured using an IniVation DVXplorer Micro event camera mounted on a stationary chaser spacecraft platform which observes a moving target spacecraft. The observed target is covered in multi-layer insulation and is equipped with a solar panel and a docking cone. The dataset contains 15 unique experiments, each approximately 300 seconds long, featuring distinct trajectories under different lighting conditions. Time-synchronized ground-truth position and velocity data are provided via a PhaseSpace motion capture system with sub-millimeter accuracy.

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

$$\mathcal{T}(x, y, t) = \exp\!\left(-\frac{t - t_{\text{last}}(x,y)}{\tau}\right)$$

where $\tau$ is a fixed time constant controlling the decay rate. This formulation is described in Lagorce et al. (2017), *"HOTS: A Hierarchy of Event-Based Time-Surfaces for Pattern Recognition,"* IEEE TPAMI. IEI-ATS replaces the fixed $\tau$ with a **per-pixel adaptive time constant** derived from local IEI statistics, and replaces the hard reset at each event with an exponential moving average (EMA) update rule.

---

## Algorithm Architecture

The processing pipeline operates on fixed-duration time windows (frames) accumulated from the event stream. Each frame passes through the following stages:

```
Events (x, y, t, p)
        │
        ▼
┌─────────────────┐
│  1. Statistics   │  Per-pixel IEI mean, std, diff
│     (+stats/)    │  Persistent EMA IEI map
└────────┬────────┘
         ▼
┌─────────────────┐
│  2. Coherence    │  Rule 1: Spatial density
│     Filtering    │  Rule 2: Temporal persistence
│  (+coherence/)   │  Rule 3: IEI regularity (CV)
│                  │  Auto-threshold (Kneedle elbow)
└────────┬────────┘
         ▼
┌─────────────────┐
│  3. Adaptive     │  Direct IEI-to-τ mapping
│     EMA          │  Attack-release envelope
│  Accumulation    │  First-order IIR surface update
│  (+accumulator/) │
└────────┬────────┘
         ▼
┌─────────────────┐
│  4. Divisive     │  Carandini-Heeger normalization
│     Normalization│  Auto-scaling semi-saturation
└────────┬────────┘
         ▼
┌─────────────────┐
│  5. Outlier      │  Per-polarity 4σ clipping
│     Rejection    │  to polarity median
└────────┬────────┘
         ▼
┌─────────────────┐
│  6. Tone Mapping │  Symmetric tanh compression
│  (+process/)     │
└────────┬────────┘
         ▼
   Output Surface
```

---

## Processing Stages

### 1. Per-Pixel IEI Statistics

Within each time window, per-pixel inter-event-interval statistics (mean, standard deviation, mean difference) are computed in batch via `computeNeighborhoodStats`. A persistent IEI map is maintained across frames using an EMA:

$$\widehat{\text{IEI}}_t(x,y) = (1 - \alpha) \cdot \widehat{\text{IEI}}_{t-1}(x,y) + \alpha \cdot \text{IEI}_t(x,y)$$

where $\alpha$ is the IEI smoothing factor (`iei_alpha`). This gives each pixel "memory" of its typical firing rate — a pixel that was active in recent frames but quiet in the current frame retains its last blended IEI value rather than snapping to zero.

### 2. Three-Rule Coherence Filtering

Three independent scoring rules produce per-pixel maps in $[0, 1]$, which are summed to form a combined coherence score. The combined score is then auto-thresholded to produce the filter mask.

**Rule 1 — Spatial density.** A KD-tree radius search in normalized $(x, y, t)$ space computes a density proxy from the summed nearest-neighbour distances within radius $r_s$. Results are log-normalized to compress the heavy-tailed distribution.

**Rule 2 — Temporal persistence.** Cross-frame persistence is evaluated by finding the nearest active pixel in the *previous* frame's trace map via a 3D KD-tree search over normalized (row, column, value) space. Pixels whose nearest predecessor exceeds `persistence_threshold` are rejected. This rule is skipped on the first frame.

**Rule 3 — IEI regularity.** The coefficient of variation (CV = $\sigma / \mu$) of the per-pixel IEI is computed directly from the existing per-pixel statistical maps (`std_map` and `mean_map`) produced in Stage 1. The regularity score is $s_{\text{reg}} = 1 / (1 + \text{CV})$. An IEI magnitude penalty prevents noise regions (where all pixels fire slowly but uniformly) from producing artificially low CV:

$$s_{\text{mag}} = \exp\!\left(-\frac{\text{IEI}}{2 \cdot \widetilde{\text{IEI}}}\right)$$

where $\widetilde{\text{IEI}}$ is the median observed IEI. The combined similarity score is $s = s_{\text{reg}} \cdot s_{\text{mag}}$. Note that this is distinct from the persistent EMA IEI map (`iei_map`) used in Stage 3: Rule 3 uses the raw within-frame IEI estimate to evaluate instantaneous spatial regularity, while the accumulator uses the smoothed multi-frame history to set decay rates.

**Auto-thresholding.** Rather than using a fixed `coherence_threshold`, the filter mask threshold is determined automatically each frame using the Kneedle algorithm (Satopaa et al., 2011, *"Finding a 'Kneedle' in a Haystack: Detecting Knee Points in System Behavior,"* Proc. ICDCSW). The blurred coherence map values are sorted to form a monotonically decreasing curve, both axes are normalized to $[0, 1]$, and the point of maximum perpendicular distance from the chord connecting the curve endpoints identifies the elbow — the natural transition between signal and noise. An EMA smooths the detected threshold across frames to prevent frame-to-frame jitter:

$$\hat{\theta}_t = \beta \cdot \theta_t + (1 - \beta) \cdot \hat{\theta}_{t-1}$$

where $\theta_t$ is the raw elbow threshold and $\beta$ is the `threshold_smoothing` parameter.

### 3. Adaptive EMA Accumulation

**Direct IEI-to-$\tau$ mapping.** The per-pixel decay time constant is set directly from the persistent IEI map:

$$\tau_{\text{active}}(x,y) = \max\!\left(\widehat{\text{IEI}}(x,y),\, \varepsilon\right)$$

followed by Gaussian spatial smoothing. A pixel firing every 0.05 s gets $\tau \approx 0.05$ s; a pixel firing every 0.5 s gets $\tau \approx 0.5$ s. This physics-based relationship is deterministic across frames — the same physical event rate always produces the same decay behaviour, regardless of what other pixels are doing.

**Attack-release envelope.** A binary activity indicator is constructed from the coherence-filtered polarity input and expanded via morphological dilation with a disk structuring element of radius 1. This dilation bridges sub-pixel gaps in sparse event data without the boundary-smearing artifacts that Gaussian smoothing of a binary mask would cause. Active pixels use $\tau_{\text{active}}$; idle pixels use the fixed release constant $\tau_{\text{release}} \ll \tau_{\text{active}}$:

$$\tau_{\text{eff}}(x,y) = \tau_{\text{active}}(x,y) \cdot m(x,y) + \tau_{\text{release}} \cdot (1 - m(x,y))$$

The activity indicator is derived from the *coherence-filtered* input (`polarity_map .* filter_mask`), not the raw polarity map. Without this ordering, pixels rejected by coherence filtering appear "active" to the envelope and receive the slow $\tau_{\text{active}}$ decay, causing trailing artifacts.

**EMA update (first-order IIR).** The blending coefficient is computed as:

$$\alpha_{\text{eff}}(x,y) = 1 - \exp\!\left(-\frac{\Delta t}{\tau_{\text{eff}}(x,y)}\right)$$

This coupling of input gain with decay rate prevents runaway accumulation. The surface update is:

$$S_t(x,y) = \alpha_{\text{eff}}(x,y) \cdot u_t(x,y) + \left(1 - \alpha_{\text{eff}}(x,y)\right) \cdot S_{t-1}(x,y)$$

where $u_t$ is the signed polarity input weighted by the continuous coherence score. After the EMA update, pixels outside the dilated activity mask are explicitly zeroed to prevent residual energy from persisting in inactive regions.

### 4. Divisive Normalization

The raw EMA surface is normalized using a variant of the canonical divisive normalization model (Carandini & Heeger, 2012, *"Normalization as a canonical neural computation,"* Nature Reviews Neuroscience, 13, 51–62):

$$\hat{S}(x,y) = \frac{S_{\text{raw}}(x,y)}{\sigma_t + C_{\text{smooth}}(x,y)^n}$$

where $C_{\text{smooth}}$ is a spatially pooled event count map (Gaussian-filtered), $n$ is a compression exponent (`div_norm_exp`, default $1.0$), and $\sigma_t$ is a per-frame dynamic scale computed as the median of the spatially smoothed EMA magnitude over all active pixels. This makes the normalization auto-scaling: $\sigma_t$ tracks the typical activity level of the scene each frame rather than requiring a manually tuned fixed value.

The normalization pool is built from unsigned event counts rather than the signed surface magnitude, preserving polarity contrast through normalization. Critically, the normalization is applied to a *separate copy* of the surface — the raw, unnormalized EMA state feeds back into the next frame's update. This prevents a feedback loop where the normalization amplifies residual values faster than the EMA decay can remove them.

### 5. Outlier Rejection

After divisive normalization, a per-polarity iterative $\sigma$-clipping step removes residual extreme values. Values beyond $n_\sigma$ standard deviations from the polarity-conditional mean are replaced with the polarity median rather than a hard bound to avoid introducing discontinuities. The clipping is applied in two passes for iterative narrowing. This handles hot pixels and isolated coherence filter leakage that survive normalization.

### 6. Tone Mapping

The final output is mapped to $[0, 1]$ using a symmetric hyperbolic tangent tone curve:

$$\hat{S}_{\text{out}}(x,y) = 0.5 + 0.5 \cdot \tanh\!\left(\frac{\hat{S}(x,y)}{s}\right)$$

where $s$ is a scalar scale parameter (default $s = 3$). Larger $s$ produces a softer curve with more headroom for extreme values; smaller $s$ increases contrast in quiet regions. A zero-valued surface maps exactly to $0.5$ (mid-gray), preserving the signed polarity symmetry of the surface. This approach follows the sigmoid tone compression operator described in Reinhard et al. (2002), *"Photographic Tone Reproduction for Digital Images,"* ACM SIGGRAPH.

---

## Repository Structure

```
IEI-ATS/
│
├── main.m                        
│
├── +accumulator/
│   └── adaptiveGlobalDecay.m 
│   └── adaptiveTimeSurfaceZhu.m 
│   └── localAdaptiveTimeSurface.m 
│   └── motionEncodedTimeSurface.m  
│   └── speedInvariantTimeSurface.m 
│   └── timeSurface.m 
│
├── +filters/
│   ├── computeCoherenceMask.m      
│   ├── findSpatialNeighbours.m     
│   ├── findSimilarities.m        
│   └── findPersistenceVectorized.m 
│   └── backgroundActivityFilter.m 
│   └── eventDensityFilter.m 
│   └── spatiotemporalCorrelation.m 
│   └── stccFilter.m 
│
├── +stats/
│   ├── computeNeighborhoodStats.m  
│   ├── findElbowThreshold.m       
│   └── printPercentComplete.m       
│
├── +process/
│   ├── generateMeshFromFrame.m
│   ├── sliceToValidRange.m  
│   ├── symmetricToneMappingNorm.m       
│   └── unwrapMap.m    
│   
└── +plot/                       

```

The repository also includes reference implementations of several published methods for benchmarking. Alternative accumulators in `+accumulator/` include the classical exponential time surface (Lagorce et al., 2017), speed-invariant time surface (Manderscheid et al., 2019), adaptive global decay (Nunes et al., 2023), and others. Alternative denoising filters in `+filter/` include the background activity filter (Delbruck, 2008), spatiotemporal correlation filter (Liu et al., 2015), event density filter (Feng et al., 2020), and the STCC-Filter (Li et al., 2024). All filters share a consistent interface and are selectable via the `filterSelection` variable in `main.m`.

---

## Dependencies

- **MATLAB R2021b or later** (uses `imdilate`, `imgaussfilt`, `imbilatfilt`)
- **Image Processing Toolbox** — required for morphological and filter operations
- **Statistics and Machine Learning Toolbox** — required for `createns`, `knnsearch`, `rangesearch`
- Event data in **HDF5 format** with datasets: `/timestamp`, `/x`, `/y`, `/polarity`

To convert AEDAT4 recordings to HDF5, see the companion import script from the [NEXUS](https://github.com/Carleton-SRL/NEXUS/tree/main/import) repository.

---

## Getting Started

### 1. Configure the data path and dataset

In `main.m`, set the path and filename for your HDF5 event recording:

```matlab
hdf5Path = '/path/to/your/datasets/';
fileName = 'your_recording.hdf5';
```

### 2. Select processing algorithms

Choose the filtering and accumulation methods:

```matlab
filterSelection      = 'COHERENCE';   % Options: 'COHERENCE', 'STC', 'BAF', 'EDF', 'STCC', 'NONE'
accumulatorSelection = 'IEI-ATS';     % Options: 'IEI-ATS', 'HOTS', 'SITS', 'METS', 'AGD', 'EVO-ATS'
```

### 3. Set the processing window

```matlab
t_start_process = 0;     % [seconds]
t_end_process   = 1000;  % [seconds]
t_interval      = 0.033; % Accumulation window [seconds]
```

### 4. Run

```matlab
main
```

Output frames are written to `videoWriters{1}` (AVI) and optionally as individual PNGs to the configured `frameOutputFolder`.

---

## Parameters

### Coherence filter

| Parameter | Description | Typical Range |
|---|---|---|
| `coh_params.r_s` | Spatial density search radius (normalized) | `0.02 – 0.08` |
| `threshold_smoothing` | EMA smoothing factor for auto-threshold | `0.0 – 1.0` |

### Accumulator

| Parameter | Description | Typical Range |
|---|---|---|
| `alts_params.dt` | Frame duration; numerator in the EMA blending coefficient | `t_interval` |
| `alts_params.filter_size` | Gaussian filter size for tau and gain smoothing | `9 – 11` |
| `alts_params.filter_sigma` | Gaussian filter sigma for tau and gain smoothing | `7.0 – 9.0` |
| `alts_params.surface_tau_release` | Release time constant for inactive pixels [s] | `1.0 – 5.0` |
| `alts_params.div_norm_exp` | Exponent $n$ controlling compression of high-activity regions | `0.5 – 1.5` |
| `iei_alpha` | IEI map EMA smoothing factor (higher = faster adaptation) | `0.5 – 0.95` |

---

## Reference Implementations

IEI-ATS includes reference implementations of several published methods for benchmarking:

**Accumulators (`+accumulator/`):**
- **Classical Time Surface** — Lagorce et al. (2017), *"HOTS: A Hierarchy of Event-Based Time-Surfaces for Pattern Recognition,"* IEEE TPAMI, 39(7), 1346–1359. [IEEE](https://doi.org/10.1109/TPAMI.2016.2574707)
- **Speed-Invariant Time Surface** — Manderscheid et al. (2019), *"Speed Invariant Time Surface for Learning to Detect Corner Points with Event-Based Cameras,"* Proc. IEEE CVPR. [IEEE](https://doi.org/10.1109/CVPR.2019.01049)
- **Adaptive Global Decay (AGD)** — Nunes et al. (2023), *"Adaptive Global Decay Process for Event Cameras,"* Proc. IEEE CVPR. [IEEE](https://doi.org/10.1109/CVPR52729.2023.00942)

**Denoising Filters (`+filter/`):**
- **Background Activity Filter (BAF)** — Delbruck (2008), *"Frame-free dynamic digital vision,"* Proc. Intl. Symp. on Secure-Life Electronics.
- **Spatiotemporal Correlation (STC)** — Liu et al. (2015), *"Design of a Spatiotemporal Correlation Filter for Event-based Sensors,"* Proc. IEEE ISCAS. [IEEE](https://doi.org/10.1109/ISCAS.2015.7168735)
- **Event Density Filter (EDF)** — Feng et al. (2020), *"Event Density Based Denoising Method for Dynamic Vision Sensor,"* Applied Sciences. [MDPI](https://doi.org/10.3390/app10062024)
- **STCC-Filter** — Li et al. (2024), *"STCC-Filter: Space-Time-Content Correlation Filter for Event Camera Noise Removal,"* Signal Processing: Image Communication. [Elsevier](https://doi.org/10.1016/j.image.2024.117135)

---

## License

MIT License. Copyright (c) 2026 Alexander Crain. See [LICENSE](LICENSE) for details.

---

## Citation

If you use IEI-ATS in your research, please cite:

```bibtex
@inproceedings{crain2026ieats,
  author    = {Crain, Alexander and Ulrich, Steve},
  title     = {Event-Based Spacecraft Representation Using Inter-Event-Interval Adaptive Time Surfaces},
  booktitle = {36th AIAA/AAS Space Flight Mechanics Meeting},
  address   = {Orlando, FL},
  month     = jan,
  year      = {2026},
  note      = {AIAA/AAS SFM 2026}
}
```

If you use the EVOS dataset, please cite:

```bibtex
@misc{crain2025evos,
  author    = {Crain, Alexander and Ulrich, Steve},
  title     = {{EVent-based Observation of Spacecraft (EVOS)}: A Neuromorphic Dataset for Spacecraft Proximity Operations},
  year      = {2025},
  publisher = {Federated Research Data Repository},
  doi       = {10.20383/103.01538},
  url       = {https://doi.org/10.20383/103.01538},
  license   = {CC BY 4.0}
}
```