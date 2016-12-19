# -*- coding: utf-8 -*-
"""
Configuration Settings
----------------------
Built for the VIRUS instrument as well as LRS2 on HET

@author: gregz
"""

Amps = ["LL", "RU"]
Amp_dict = {"LL": ["LU","L"], "RU": ["RL","R"]}
virus_wl = {"LL": [3490,5500], "RU": [3490,5500]}
virus_sn = {"LL": "virus", "RU": "virus"}
lrs2b_wl = {"LL": [3633,4655], "RU": [4550,7000]}
lrs2b_sn = {"LL": "lrs2_uv", "RU": "lrs2_orange"}
lrs2r_wl = {"LL": [3490,5500], "RU": [3490,5500]}
lrs2r_sn = {"LL": "lrs2_red", "RU": "lrs2_farred"}
virus_fs = 8.
lrs2b_fs = 6.
lrs2r_fs = 6.
virus_di = {"LL": 1.9, "RU": 1.9}
lrs2b_di = {"LL": 0.5, "RU": 1.2}
lrs2r_di = {"LL": 1.2, "RU": 1.2}