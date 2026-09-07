import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import muspan as ms
import shift

#%% Generate example datasets
# MuSpAn domains with A-B local correlation, C random
# MuSpAn domains with A-B local exclusion, C random
np.random.seed(0)
domains = []
n_points_per_population = 500
n_seeds = 20
for i in range(10):
    d = ms.domain(f"State1_{i}")
    print(d.name)

    seeds = np.random.uniform(0, 1000, (n_seeds, 2))

    A = np.clip(seeds[np.random.choice(n_seeds, n_points_per_population)] + np.random.normal(0, 30, (n_points_per_population, 2)), 0, 1000)
    B = np.clip(seeds[np.random.choice(n_seeds, n_points_per_population)] + np.random.normal(0, 30, (n_points_per_population, 2)), 0, 1000)
    C = np.random.uniform(0, 1000, (n_points_per_population, 2))

    d.add_points(np.vstack([A, B, C]))
    d.add_labels("Celltype", ["A"] * n_points_per_population + ["B"] * n_points_per_population + ["C"] * n_points_per_population)
    bdy = np.array([[0,0],[0,1000],[1000,1000],[1000,0]])
    d.estimate_boundary(method='specify',specify_boundary_coords=bdy)
    domains.append(d)

for i in range(10):
    d = ms.domain(f"State2_{i}")
    print(d.name)

    seeds = np.random.uniform(0, 1000, (n_seeds, 2))

    A = np.clip(seeds[np.random.choice(n_seeds, n_points_per_population)] + np.random.normal(0, 30, (n_points_per_population, 2)), 0, 1000)

    B = np.empty((0, 2))
    while len(B) < n_points_per_population:
        pts = np.random.uniform(0, 1000, (2*n_points_per_population, 2))
        keep = np.min(np.linalg.norm(pts[:, None, :] - seeds[None, :, :], axis=2), axis=1) > 100
        B = np.vstack([B, pts[keep]])
    B = B[:n_points_per_population]

    C = np.random.uniform(0, 1000, (n_points_per_population, 2))

    d.add_points(np.vstack([A, B, C]))
    d.add_labels("Celltype", ["A"] * n_points_per_population + ["B"] * n_points_per_population + ["C"] * n_points_per_population)

    bdy = np.array([[0,0],[0,1000],[1000,1000],[1000,0]])
    d.estimate_boundary(method='specify',specify_boundary_coords=bdy)
    domains.append(d)
    
    
    
#%% CALCULATE CONSISTENT STATISTICS FOR EACH DOMAIN, COMBINE INTO DATAFRAME
# statistics = ["counts", "persistent homology", "quadrat correlation matrix", "wasserstein distance", "cross PCF", "topographical correlation map level set filtration", "adjacency permutation test"]
# args_for_statistics = 
rows = {}
for domain in domains:
    print(domain.name)
    df_row = shift.generate_feature_dataframe(domain, 'Celltype')
    df_row["name"] = domain.name
    df_row["State"] = domain.name.split("_")[0]
    rows[domain.name] = df_row
df = pd.concat(rows, ignore_index=True)
df = df[["name", "State"] + [c for c in df.columns if c not in ("name", "State")]]

# Optional - save/load if you want to calculate SHIFT in a subsequent workflow
# df.to_csv('./Example_Statistics_Dataframe.csv', index=False)
# df = pd.read_csv('./Example_Statistics_Dataframe.csv')
#%% CALCULATE SHIFT SCORE
SHIFT, pvals = shift.compute_shift(df, comparator_state_column="State", metadata_columns=["name"], return_component_pvalues=True)




