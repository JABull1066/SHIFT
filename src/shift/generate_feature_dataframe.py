import numpy as np
import muspan as ms
import pandas as pd

DEFAULT_STATISTICS = [
    "counts",
    "quadrat correlation matrix",
    "wasserstein distance",
    "cross PCF"
]

def generate_feature_dataframe(domain, cell_label_name,  statistics=None, min_cells_required=20, args_for_statistics=None):
    """
    Calculate statistics for a given domain and set of cell pairs.

    Parameters
    ----------
    domain : muspan.Domain
        The spatial domain containing the cells.

    cell_label_name : str
        The name of a categorical MuSpAn label within the domain, used to determine which analyses to include.

    statistics : list[str]
        Which statistics to calculate. Supported values include:
        - "counts": cell counts of each cell type
        - "persistent homology": mean interaction strength between cell pairs.
        If None, defaults to all supported statistics.

    min_cells_required : int
        Minimum number of cells of a given type required to calculate statistics for that type. Default is 20.

    args_for_statistics : dict, optional
        Additional arguments for specific statistics. For example, for "persistent homology", this could include 'threshold_radius' to specify the radius of a loop for nLoops in the Vietoris-Rips filtration.

    Returns
    -------
    dict[str, float]
        A dictionary mapping statistic names to their calculated values.
    """

    results = {}

    if cell_label_name not in domain.labels:
        raise ValueError(f"Label '{cell_label_name}' not found in domain labels.")

    if statistics is None:
        statistics = DEFAULT_STATISTICS

    if args_for_statistics is None:
        args_for_statistics = {
            "persistent_homology": {},
            "quadrat_correlation_matrix": {},
            "cross PCF": {},
            "topographical_correlation_map": {},
            "adjacency_permutation_test": {},
            "morisita_horn_index": {}
        }


    cell_pop = ms.query.query(domain, ('label',cell_label_name),'is not',None)
    counts, cts = ms.summary_statistics.label_counts(domain, cell_label_name, normalised=False)
    if "counts" in statistics:
        for i, celltype in enumerate(cts):
            results[f"{celltype}-None_Summaries_Count"] = counts[i]

    if "persistent homology" in statistics:
        for i, celltype in enumerate(cts):
            if float(counts[i]) > min_cells_required:
                # Standard PH on this point population alone
                TDA_dict = ms.topology.vietoris_rips_filtration(domain,population=(cell_label_name,celltype))
                vec, stats = ms.topology.vectorise_persistence(TDA_dict)
                for j in range(len(vec)):
                    statname = '-'.join(stats[j].split(' '))
                    results[f'{celltype}-None_PH_{statname}'] = vec[j]
                results[f'{celltype}-None_PH_nLoops'] = _get_n_loops(TDA_dict['dgms'][1], args_for_statistics.get('persistent_homology', {}).get('threshold_radius', 20))


    if "quadrat correlation matrix" in statistics:
        side_length = args_for_statistics.get('quadrat_correlation_matrix', {}).get('side_length', 100)
        SES, A, cats = ms.region_based.quadrat_correlation_matrix(domain,cell_label_name,population=cell_pop, region_kwargs={'side_length':side_length},low_observation_bound=min_cells_required)
        for i, cat_i in enumerate(cats):
            for j, cat_j in enumerate(cats):
                if i != j:
                    results[f'{cat_i}-{cat_j}_QCM_SES'] = SES[i,j]

    if "wasserstein distance" in statistics:
        for i, cat_i in enumerate(cats):
            for j, cat_j in enumerate(cats):
                if i != j:
                    w = ms.distribution.sliced_wasserstein_distance(domain, (cell_label_name, cat_i), (cell_label_name, cat_j))
                    results[f'{cat_i}-{cat_j}_Wasserstein_WassersteinDistance'] = w

    if "cross PCF" in statistics:
        for i, cat_i in enumerate(cats):
            for j, cat_j in enumerate(cats):
                if float(counts[i]) > min_cells_required and float(counts[j]) > min_cells_required:
                    max_R = args_for_statistics.get('cross PCF', {}).get('max_R', 150)
                    annulus_step = args_for_statistics.get('cross PCF', {}).get('annulus_step', 5)
                    annulus_width = args_for_statistics.get('cross PCF', {}).get('annulus_width', 10)
                    r, g = ms.spatial_statistics.cross_pair_correlation_function(domain, (cell_label_name, cat_i), (cell_label_name, cat_j), max_R=max_R, annulus_step=annulus_step, annulus_width=annulus_width, visualise_output=False)
                    results[f'{cat_i}-{cat_j}_CrossPCF_gmax'] = g[np.nanargmax(g)]
                    results[f'{cat_i}-{cat_j}_CrossPCF_rpeak'] = r[np.nanargmax(g)]
                    results[f'{cat_i}-{cat_j}_CrossPCF_gmin'] = g[np.nanargmin(g)]
                    results[f'{cat_i}-{cat_j}_CrossPCF_rtrough'] = r[np.nanargmin(g)]
                    for mult in [1, 3, 5, 10, 20]:
                        # Take the 1st, 3rd, 5th, 10th, and 20th multiples of the annulus width for integration
                        rmax = mult * annulus_width
                        if rmax > max_R:
                            continue
                        idx = np.searchsorted(r, rmax, side="right")
                        results[f'{cat_i}-{cat_j}_CrossPCF_int-0-{rmax}'] = np.trapezoid(g[:idx], r[:idx])

    if "topographical correlation map level set filtration" in statistics:
        for i, cat_i in enumerate(cts):
            pop_a = ms.query.query(domain, ('label',cell_label_name),'is',cat_i)
            if float(counts[i]) <= min_cells_required:
                continue
            for j, cat_j in enumerate(cts):
                pop_b = ms.query.query(domain, ('label',cell_label_name),'is',cat_j)
                if float(counts[j]) <= min_cells_required:
                    continue
                radius_of_interest = args_for_statistics.get('topographical_correlation_map', {}).get('radius_of_interest', 50)
                kernel_radius = args_for_statistics.get('topographical_correlation_map', {}).get('kernel_radius', 150)
                kernel_sigma = args_for_statistics.get('topographical_correlation_map', {}).get('kernel_sigma', 50)
                TCM = ms.spatial_statistics.topographical_correlation_map(domain, pop_a, pop_b,radius_of_interest=radius_of_interest, kernel_radius=kernel_radius, kernel_sigma=kernel_sigma,visualise_output=False)
                lsf = ms.topology.level_set_filtration(TCM,visualise_output=False)
                vec, stats = ms.topology.vectorise_persistence(lsf)
                for j in range(len(vec)):
                    statname = '-'.join(stats[j].split(' '))
                    results[f'{cat_i}-{cat_j}_TCM-LS_{statname}'] = vec[j]

    if "adjacency permutation test" in statistics:

        ms.networks.generate_network(domain, 'network_name_temp', objects_as_nodes=cell_pop, **args_for_statistics.get('adjacency_permutation_test', {}).get('generate_network', {}))
        SES, A, cats = ms.networks.adjacency_permutation_test(domain, 'network_name_temp', cell_label_name, population=cell_pop)
        for i, cat_i in enumerate(cats):
            for j, cat_j in enumerate(cats):
                results[f'{cat_i}-{cat_j}_APT_SES'] = SES

    if "morisita horn index" in statistics:
        SES, A, cats = ms.region_based.morisita_horn_index(domain, cell_label_name, region_kwargs={'side_length':200})
        for i, cat_i in enumerate(cats):
            for j, cat_j in enumerate(cats):
                if i != j:
                    results[f'{cat_i}-{cat_j}_MorisitaHorn_SES'] = SES[i,j]

    return pd.DataFrame([results])


def _get_n_loops(dgm, threshold_radius=20):
    rs = dgm[:,1]-dgm[:,0]
    
    # Critical r
    critical_r = threshold_radius*np.sqrt(3)
    plotmask = rs > critical_r
    n_loops = np.sum(plotmask)

    return n_loops