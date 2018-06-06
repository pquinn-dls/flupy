├── Acquisition_instrument
│   ├── SEM
│   │   ├── Detector
│   │   │   └── EDS
│   │   │       ├── azimuth_angle (º)
│   │   │       ├── elevation_angle (º)
│   │   │       ├── energy_resolution_MnKa (eV)
│   │   │       ├── live_time (s)
│   │   │       └── real_time (s)
│   │   ├── beam_current (nA)
│   │   ├── beam_energy (keV)
│   │   ├── convergence_angle (mrad)
│   │   ├── magnification
│   │   ├── microscope
│   │   ├── Stage
│   │   │   ├── rotation (º)
│   │   │   ├── tilt_alpha (º)
│   │   │   ├── tilt_beta (º)
│   │   │   ├── x (mm)
│   │   │   ├── y (mm)
│   │   │   └── z (mm)
│   │   └── working_distance (mm)
│   └── TEM
│   │   ├── Detector
│   │   ├   ├── EDS
│   │   ├   │   ├── azimuth_angle (º)
│   │   │   │   ├── elevation_angle (º)
│   │   │   │   ├── energy_resolution_MnKa (eV)
│   │   │   │   ├── live_time (s)
│   │   │   │   └── real_time (s)
│   │   │   └── EELS
│   │   │       ├── aperture (mm)
│   │   │       ├── collection_angle (mrad)
│   │   │       ├── dwell_time (s)
│   │   │       ├── exposure (s)
│   │   │       ├── frame_number
│   │   │       └── spectrometer
│   │   ├── Biprism
│   │   │   ├── azimuth_angle (º)
│   │   │   ├── position
│   │   │   └── voltage (V)
│   │   ├── acquisition_mode
│   │   ├── beam_current (nA)
│   │   ├── beam_energy (keV)
│   │   ├── camera_length (mm)
│   │   ├── convergence_angle (mrad)
│   │   ├── magnification
│   │   ├── microscope
│   │   └── Stage
│   │       ├── rotation (º)
│   │       ├── tilt_alpha (º)
│   │       ├── tilt_beta (º)
│   │       ├── x (mm)
│   │       ├── y (mm)
│   │       └── z (mm)
│   ├── XRF
│   │   ├── incident_angle (º)
│   │   ├── exit_angle (º)
│   │   ├── Detector
│   │   │   ├── energy_resolution_MnKa (keV)
│       │   │   ├── gain (keV/ch)
│       │   │   ├── offset (keV)
│       │   │   ├── fano 
│       │   │   ├── live_time (s)
│       │   │   ├── real_time (s)
│       │   │   ├── detector_distance (mm)
│       │   │   ├── detector_area (mm)
│       │   │   ├── detector_type (Si or Ge)
│       │   │   ├── no_detector_sensors
│       │   │   ├── attenuators 
│       │   │   │   ├── window
│       │   │   │   │  ├─ thickness
│       │   │   │   │  ├─ composition
│       │   │   │   ├── sensor    
│       │   │   │   │  ├─ thickness
│       │   │   │   │  ├─ composition
│       │   │   │   ├── deadlayer
│       │   │   │   │  ├─ thickness
│       │   │   │   │  ├─ composition
│       │   │   │   ├── contacts
│       │   │   │   │  ├─ thickness
│       │   │   │   │  ├─ composition
│       │   ├── Attenuators      *thickness, composition pairs - any arbitrary description key
│       │   │   ├── freespace
│       │   │   │  ├─ thickness
│       │   │   │  ├─ composition
│       │   │   ├── sample_window
│       │   │   │  ├─ thickness
│       │   │   │  ├─ composition
│       │   └── XRFCalibrationStandard
│       │   ├── beam_energy (keV)
│       ├── XRD **** Future ****
│       └── Stage
│           ├── rotation (º)
│           ├── tilt_alpha (º)
│           ├── tilt_beta (º)
│           ├── x (mm)
│           ├── y (mm)
│           └── z (mm)
├── General
│   ├── authors
│   ├── date
│   ├── doi
│   ├── original_filename
│   ├── notes
│   ├── time
│   ├── time_zone
│   └── title
├── Sample
│   ├── credits
│   ├── description
│   ├── elements
│   ├── xray_lines
│   ├── composition
│   ├── concentrations        
│   ├── matrix        
│   └── multilayer     *For XRF calculations with thick samples - absorption effects *
└── Signal
    ├── Noise_properties
    │   ├── Variance_linear_model
    │   │   ├── correlation_factor
    │   │   ├── gain_factor
    │   │   ├── gain_offset
    │   │   └── parameters_estimation_method
    │   └── variance
    ├── binned
    ├── quantity
    ├── signal_type
    └── signal_origin