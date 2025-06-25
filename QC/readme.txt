Beam angle optimization for three cases:
1. abdomen (1243050)
2. lung (7119049)
3. HN02

For each case, we have 72 available angles: (24 gantry angles (0, 15, ..., 345) x 3 couch angles (0, 30, 60).

Out of 72 angles, we need to choose (i) 4 angles for HN02 and (ii) 3 angles for abdomen and lung case.

Following files are included in the shared folder:
1. {case id}_{gantry angle}_{couch angle}.mat: Dij information for the corresponding gantry and couch angles
2. {case id}.mat: file containing cst, ct info
3. "test_admm_..." files: these files contain information about dvh constraints, objective l2 term weights, hyperparameter values, etc. for the test cases.
