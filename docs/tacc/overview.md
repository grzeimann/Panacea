# TACC Overview

This page is for HET users running or accessing Panacea reductions on the Texas Advanced Computing Center (TACC).

## Accounts and Access
- Create an account: https://portal.tacc.utexas.edu/
- After your account is created, email Greg Zeimann <gregz@astro.as.utexas.edu> your TACC username to be added to the HET group.
- SSH into TACC:
```bash
ssh -Y USERNAME@ls6.tacc.utexas.edu
```

## Where to find data products
Daily automated reductions are written to:
```text
/work/03946/hetdex/maverick/LRS2/PROGRAM-ID
```
Replace PROGRAM-ID with your program number (e.g., HET19-1-999).

To copy reductions for your program to your local machine:
```bash
scp -r username@ls6.tacc.utexas.edu:/work/03946/hetdex/maverick/LRS2/PROGRAM-ID .
```

To fetch specific files for a given date:
```bash
scp username@ls6.tacc.utexas.edu:/work/03946/hetdex/maverick/LRS2/PROGRAM-ID/spec*20190105*.fits .
```

See also: [Running on TACC](./running.md), [Data Products](../data-products/overview.md)
