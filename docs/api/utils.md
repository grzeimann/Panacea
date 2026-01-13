# panacea.utils

This page documents selected public utilities in panacea.utils that are commonly used in calibration and extraction workflows.

```{eval-rst}
.. automodule:: panacea.utils
   :members: safe_division, build_weight_matrix, robust_polyfit, read_arc_lines, get_config_file, get_bigarray, find_lines, create_header_objection
   :undoc-members:
   :show-inheritance:
```

Examples
- Safe division guarding against zeros/NaNs:

```python
import numpy as np
from panacea.utils import safe_division

num = np.array([1.0, 2.0, -3.0, 4.0])
denom = np.array([1.0, 0.0, np.inf, 1e-12])
out = safe_division(num, denom, eps=1e-8, fillval=0.0)
print(out)  # array([1.0, 0.0, 0.0, 0.0])
```

- Build a spatial weight matrix between fibers:

```python
import numpy as np
from panacea.utils import build_weight_matrix

x = np.array([0.0, 1.0, 0.0, 1.0])
y = np.array([0.0, 0.0, 1.0, 1.0])
W = build_weight_matrix(x, y, sig=0.75)
# columns sum to 1.0
assert np.allclose(W.sum(axis=0), 1.0)
```
