# Tumor Invasion Model: Leader–Follower Dynamics using CompuCell3D

This repository contains simulation code and data for modeling biophysical regulation of leader–follower heterogeneity in collective tumor invasion, using the Cellular Potts Model framework implemented in CompuCell3D.

---

## 📦 Installation

This model runs on **CompuCell3D (v4.3.1)**.

### ✅ 1. Download & Install CompuCell3D

Visit the [CompuCell3D download page](https://compucell3d.org/Downloads) and install for your OS:

- Windows: Installer available
- macOS/Linux: Use conda installation

Or install using `conda`:

```bash
conda create -n cc3d-env -c compucell3d -c conda-forge compucell3d=4.3.1
conda activate cc3d-env
```

### ✅ 2. Clone this Repository

```bash
git clone https://github.com/your-username/tumor-invasion-model.git
cd tumor-invasion-model
```

---

## 🧪 Running the Simulation

The model uses multiple `Steppables` defined in `CCIecmSteppables.py` and can be launched using:

```python
from cc3d import CompuCellSetup
from CCIecmSteppables import (
    NeighborTrackerPrinterSteppable,
    ConstraintInitializerSteppable,
    GrowthSteppable,
    MitosisSteppable
)

CompuCellSetup.register_steppable(NeighborTrackerPrinterSteppable(frequency=100))
CompuCellSetup.register_steppable(ConstraintInitializerSteppable(frequency=1))
CompuCellSetup.register_steppable(GrowthSteppable(frequency=1))
CompuCellSetup.register_steppable(MitosisSteppable(frequency=1))

CompuCellSetup.run()
```

Then run:

```bash
cc3d-run -i tumor_invasion_project.cc3d
```

---

## 🔁 Parameter Space

```python
parameter_values = {
    'Jlf': [-5 to 5],
    'mu': [0 to 30],
    'PP': [0.0 to 1.0]
}
```

Each simulation is indexed by `iteration`, used to load the corresponding parameter combination.

---

## 📁 Output Includes

- `Metrics_Data_*.csv`: Invasive/finger area, clusters, defectors
- `ClusterComposition_*.csv`: Composition of each detached cluster
- `.png`: Boundary and phenotype plots
- Subfolder for tumor position, leader/follower positions

---

## 📜 License

MIT License

---

## 📬 Contact

Email: sheriffakeeb@gmail.com
