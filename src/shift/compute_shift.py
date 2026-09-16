import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, chi2

def compute_shift(
    df,
    comparator_state_column,
    descriptor_columns=None,
    metadata_columns=None,
    celltypes_of_interest=None,
    return_component_pvalues=False,
    states_to_compare=None
):
    """
    Compute SHIFT score from ROI-level descriptors.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per ROI. Must contain descriptor columns and a comparator state column.

    comparator_state_column : str
        Column name of the comparator state. Should contain two unique values, e.g. "control" and "treatment".

    descriptor_columns : list[str], optional
        Which columns to use for SHIFT calculation. Note all descriptor columns should be labelled "X-Y_StatisticClass_StatisticName" (e.g. "Tcell-Macrophage_CrossPCF_rpeak"). For single cell statistics, the second celltype should be "None" or the same as the first celltype (e.g. "Tcell-None_Summaries_Count").
        If None, defaults to all columns (excluding the comparator state column).

    metadata_columns : list[str], optional
        Which columns to treat as metadata (not used in SHIFT calculation).

    celltypes_of_interest : list[str], optional
        If provided, only descriptors involving these cell types will be used for SHIFT calculation. If None, defaults to all cell types present in the descriptor columns.

    return_component_pvalues : bool, optional
        If True, return the components of the SHIFT score - i.e., p-values for each `descriptor_column` and the overall SHIFT score. Default is False.

    states_to_compare : list[str], optional
        If provided, only ROIs with these two states will be used for SHIFT calculation. Must contain exactly two strings.

    Returns
    -------
    pandas.Series or numpy.ndarray
        SHIFT score per ROI.
    """

    if descriptor_columns is None:
        descriptor_columns = [col for col in df.columns if col != comparator_state_column and col not in (metadata_columns or [])]

    if states_to_compare is not None:
        if len(states_to_compare) != 2:
            raise ValueError("`states_to_compare` must contain exactly two unique values.")
        df = df[df[comparator_state_column].isin(states_to_compare)]

    X = df[descriptor_columns]

    if celltypes_of_interest is not None:
        # Filter descriptor columns to only include those involving the specified cell types
        filtered_columns = []
        for col in descriptor_columns:
            celltype1, celltype2 = col.split('_')[0].split('-')
            if celltype1 in celltypes_of_interest and celltype2 in celltypes_of_interest:
                filtered_columns.append(col)
        X = X[filtered_columns]
    else:
        celltypes_of_interest = np.unique([col.split('_')[0].split('-')[0] for col in descriptor_columns] + [col.split('_')[0].split('-')[1] for col in descriptor_columns])

    states = df[comparator_state_column]
    unique_states = states.dropna().unique()
    assert len(unique_states) == 2, f"Expected exactly two unique values in column '{comparator_state_column}', got {len(unique_states)}. If necessary, use the states_to_compare argument to specify which two states to compare."

    # Get all p values
    print(f"State 1: {unique_states[0]}, State 2: {unique_states[1]}")
    state_1_mask = df[comparator_state_column] == unique_states[0]
    state_2_mask = df[comparator_state_column] == unique_states[1]

    stats = X.columns
    pvals = [_getPvals(X, state_1_mask, state_2_mask, stat) for stat in stats]

    # Now we calculate the SHIFT score
    testdf = pd.DataFrame(pvals)
    # Sort dataframe so that all ct1 == ct2 are at the top; this is useful for drop_duplicates in a minute, as otherwise single-cell metrics will potentially be influenced by ROIs that have been dropped due to insufficient numbers of the second cell
    samecell_mask = testdf.ct1 == testdf.ct2
    testdf = pd.concat([testdf[samecell_mask], testdf[~samecell_mask]], axis=0).reset_index(drop=True)

    testdf['pvalue'][testdf['pvalue'].isna()] = 1
    testdf['effect size'] = testdf['U']/(testdf['n_state1_ROIs']*testdf['n_state2_ROIs'])
    testdf['rank-biserial correlation'] = 1 - 2*testdf['U']/(testdf['n_state1_ROIs']*testdf['n_state2_ROIs'])

    # # # Adjust p values
    # # from scipy.stats import false_discovery_control
    # testdf['log p'] = np.log10(testdf['pvalue'])

    # # # Use sign of rank-biserial correlation
    # testdf['log p directional'] = testdf['log p']
    # flip_signs = np.sign(testdf['rank-biserial correlation']) == 1
    # testdf['log p directional'][flip_signs] = -1*testdf['log p directional'][flip_signs]


    # Right, now calculate SHIFT for each cell pair of interest
    # First, get all ordered pairs of cell types in the list celltypes_of_interest
    pairs = [(ct1, ct2) for ct1 in celltypes_of_interest for ct2 in celltypes_of_interest]
    shift_scores = []
    for pair in pairs:
        ct1, ct2 = pair
        SHIFT = _get_shift_for_cellpair(testdf, ct1, ct2)
        shift_scores.append({'ct1':ct1, 'ct2':ct2, 'SHIFT': SHIFT})

    shift_scores = pd.DataFrame(shift_scores)
    if return_component_pvalues:
        return shift_scores, pvals
    else:
        return shift_scores


def _getPvals(X, state_1_mask, state_2_mask, stat):
    celltype1, celltype2 = stat.split('_')[0].split('-')
    state_1_vals = np.array(X[state_1_mask][stat])
    state_2_vals = np.array(X[state_2_mask][stat])
    U, p = mannwhitneyu(x=state_1_vals, y=state_2_vals, nan_policy='omit')
    dataframe_row = {'ct1':celltype1,'ct2':celltype2,'Statistic':stat,'U':U,'pvalue':p,'n_state1_ROIs':len(state_1_vals), 'n_state2_ROIs':len(state_2_vals),'median_state1':np.median(state_1_vals), 'median_state2':np.median(state_2_vals)}
    return dataframe_row




def _get_shift_for_cellpair(all_statistics_df, ct1, ct2):

    mask = ((all_statistics_df.ct1 == ct1) & (all_statistics_df.ct2 == ct2))
    if sum(mask) == 0:
        return 0
    filtered_df = all_statistics_df[mask]

    # Use Fishers Method
    pvalue_dataframe = (
        filtered_df
            .drop_duplicates(subset="Statistic", keep="first")
            .assign(ID="single")
            .pivot(index="ID", columns="Statistic", values="pvalue")
    )
    sig_score_fisher = float(-2*np.sum(np.log(pvalue_dataframe),axis=1)['single'])
    k = np.shape(pvalue_dataframe)[1]# Number of tests being combined
    pval = 1 - chi2.cdf(sig_score_fisher,2*k)
    
    # Now assign signs
    # Pull out specific statistics
    # {stat subset : direction of stat change associated with increased order}
    # Note that each of these contains an underscore, which is protected in the statistic name, so we can use this to identify the statistic type
    stats_for_sign_check = {
                        'QCM_SES' : 1,
                        'APT_SES' : 1,
                        'MorisitaHorn_SES' : 1,
                        'PCF_int' : 1,
                        'PCF_gmax' : 1, 
                        'PCF_gmin' : 1,
                        'PCF_rtrough' : 1,
                        'PCF_rpeak' : -1,
                        'Wasserstein_WassersteinDistance' : -1
                        }
    consensus = []
    stats = list(filtered_df.Statistic)
    for stat in stats:
        for sc in stats_for_sign_check:
            if sc in stat:
                x = filtered_df[filtered_df['Statistic'] == stat]
                sign = np.sign(float(x['rank-biserial correlation'].iloc[0])) == np.sign(stats_for_sign_check[sc])
                consensus.append(sign)
                continue
    
    
    if np.sum(consensus) > 0.5*len(consensus):
        sign = 1
    else:
        sign = -1
    # Cap SHIFT scores at 20
    val = max([np.log10(pval),np.float64(-20)])
    SHIFT = -sign*val
    return SHIFT