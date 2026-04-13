# COBRApy Model Support

This fork adds the ability to pass COBRApy (`cobra.Model`) objects directly to SMETANA, bypassing the need for SBML file paths.

## Installation

Requires both `smetana` and `cobra`:

```bash
pip install cobra
```

## Usage

### Pass COBRApy models to `main()`

```python
import cobra
from smetana.interface import main

m1 = cobra.io.read_sbml_model("species1.xml")
m2 = cobra.io.read_sbml_model("species2.xml")

# Pass COBRApy models directly -- conversion is automatic
main([m1, m2], mode="global", output="results",
     media="M9", mediadb="media_db.tsv")
```

### Convert explicitly

```python
from smetana import convert_cobrapy_model

reframed_model = convert_cobrapy_model(cobra_model)
```

### Convert a batch

```python
from smetana import convert_cobrapy_models

reframed_models = convert_cobrapy_models([cobra_model_1, cobra_model_2])
```

### Use `CobraModelCache` for repeated access

`CobraModelCache` converts all models once upfront and implements the same interface as reframed's `ModelCache`:

```python
from smetana import CobraModelCache

cache = CobraModelCache([cobra_model_1, cobra_model_2])
cache.get_ids()                        # ['species1', 'species2']
cache.get_model('species1')            # returns reframed CBModel
```

## What the conversion does

`convert_cobrapy_model()` handles four things that reframed's built-in `from_cobrapy()` does not:

1. **ID prefix normalization** -- Adds `R_` to reaction IDs and `M_` to metabolite IDs to match the BiGG namespace that SMETANA's internals expect.
2. **External compartment detection** -- Marks the extracellular compartment (tries `e`, `e0`, `extracellular`, then falls back to boundary-reaction frequency analysis).
3. **Reaction type classification** -- Assigns `EXCHANGE`, `SINK`, `TRANSPORT`, or `ENZYMATIC` types based on stoichiometric structure.
4. **Bound normalization** -- Converts COBRApy's literal +/-1000 bounds to +/-inf, matching reframed's SBML loader behavior.

## Batch script

`scripts/run_kchip_smetana.py` demonstrates batch usage across model permutations and JSON-defined media:

```bash
# Sequential
python scripts/run_kchip_smetana.py -v --output-dir results

# Parallel (8 workers)
python scripts/run_kchip_smetana.py -w 8 --output-dir results
```
