# panacea.astrometry

Astrometry utilities for converting between IFU slot coordinates and sky coordinates.

```{eval-rst}
.. automodule:: panacea.astrometry
   :members: Astrometry, IFU, FPlane
   :undoc-members:
   :show-inheritance:
```

Example
- Initialize an Astrometry model and query an IFU slot position:

```python
from panacea.astrometry import Astrometry

astro = Astrometry(ra0=150.0, dec0=2.2, pa=0.0, x0=0.0, y0=0.0)
ra, dec = astro.get_ifuslot_ra_dec('054')
print(ra, dec)
```
