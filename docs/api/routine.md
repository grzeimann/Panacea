# panacea.routine

High-level routines used in calibration and source extraction.

```{eval-rst}
.. automodule:: panacea.routine
   :members: fit_response_cont, get_response, extract_source, find_source
   :undoc-members:
   :show-inheritance:
```

Examples
- Estimate a smooth continuum for a noisy spectrum:

```python
import numpy as np
from panacea.routine import fit_response_cont

wv = np.linspace(4000, 5000, 501)
sky = 1.0 + 0.001*(wv - wv.mean())
sky[250] += 0.8  # narrow emission feature
cont = fit_response_cont(wv, sky, fil_len=21)
```

- Compute a response vector from a recognized standard (returns None if not found):

```python
import numpy as np
from panacea.routine import get_response

wv = np.linspace(4000, 5000, 501)
sky = 1.0 + 0.001*(wv - wv.mean())
resp = get_response('HZ_44', wv, sky, 'red')
if resp is not None:
    calibrated = sky * resp
```
