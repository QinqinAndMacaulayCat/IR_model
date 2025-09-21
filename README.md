
# IR_model

A Python package for modeling interest rates using equilibrium term structure models, supporting multiple short rate models including Vasicek, Hull-White, and CIR.


## Installation

Clone the repository:

```bash
git clone https://github.com/QinqinAndMacaulayCat/IR_model.git
cd IR_model
```

### Using Poetry (recommended)

```bash
poetry install
```

### Using pip

```bash
pip install -r requirements.txt
```

For development/editable install:

```bash
pip install -e .
```



## Features

- **Multiple Models**: Simulate, calibrate, and analyze classic short rate models (Vasicek, Hull-White, CIR).
- **Parameter Recovery**: Tools for MLE and risk-neutral parameter recovery.
- **Extensible Design**: Easily add new models by extending the `BaseModel` class.
- **Data Utilities**: Fetch and use sample datasets for model calibration and testing.



## Models Supported

- **Vasicek**
- **Hull-White**
- **Cox-Ingersoll-Ross (CIR)**

See [`doc/model.md`](doc/model.md) for mathematical details and model theory.

## Project Structure

```
irmodel/
	__init__.py
	data/
		data_fetcher.py
		*.rda
	model/
		__init__.py
		BaseModel.py
		Vasicek.py
		Hull_White.py
		CIR.py
scripts/
	Vasicek_Analysis.py
	HullWhite_Analysis.py
	CIR_analysis.py
```


## Usage

To run a model analysis script:

```bash
cd scripts
python Vasicek_Analysis.py   # or HullWhite_Analysis.py, CIR_analysis.py
```

Or from the project root:

```bash
PYTHONPATH=. python scripts/Vasicek_Analysis.py
```


## Extending

- Add new models in `irmodel/model/` by subclassing `BaseModel`.
- Add new data fetchers in `irmodel/data/`.


## References

- See `doc/model.md` for model theory and formulas.
- See `scripts/` for example analyses.