# Examples

This page collects short, copy-pasteable examples for common Panacea tasks.

Prerequisites
- Install Panacea and dev extras for the CLI and examples:

```bash
pip install .[dev]
```

CLI: Reduce a small dataset
- Run the pipeline on a directory or tarball. Use --help to see options.

```bash
panacea-lrs2 --help
panacea-lrs2 --input /path/to/lrs2/night/ --channel red --output ./out
```

- Tips
  - Use --config to point to a YAML config if needed.
  - Logs will indicate which frames are used and where products are written.

Python API: Smooth a sky spectrum continuum

```python
import numpy as np
from panacea.routine import fit_response_cont

wv = np.linspace(4000.0, 5000.0, 501)
# synthetic spectrum with a gentle slope and a narrow line
sky = 1.0 + 0.001 * (wv - wv.mean())
sky[250] += 0.8
cont = fit_response_cont(wv, sky, fil_len=21)
```

Python API: Safe numerical utilities

```python
import numpy as np
from panacea.utils import safe_division, build_weight_matrix

num = np.array([1.0, 2.0, -3.0, 4.0])
denom = np.array([1.0, 0.0, np.inf, 1e-12])
out = safe_division(num, denom, eps=1e-8, fillval=0.0)

x = np.array([0.0, 1.0, 0.0, 1.0])
y = np.array([0.0, 0.0, 1.0, 1.0])
W = build_weight_matrix(x, y, sig=0.75)
```

Python API: Astrometry quicklook

```python
from panacea.astrometry import Astrometry
astro = Astrometry(ra0=150.0, dec0=2.2, pa=0.0, x0=0.0, y0=0.0)
ra, dec = astro.get_ifuslot_ra_dec('054')
```

See also
- API reference: ../api/index.md
- CLI Usage: ./cli.md
- Configuration: ./configuration.md
