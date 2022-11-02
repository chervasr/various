def get_false_discovery_rate(pval, column):
    """
    _summary_

    Parameters
    ----------
    pval : pd.DataFrame
        _description_
    column: str | int
        _description_

    Returns
    -------
    list
        _description_
    """
    p_adj_list = []

    nvalues = len(pval)
    pval = pval.sort_values(by=column)
    for ix, p in enumerate(pval[column]):
        fdr_adj_p_val = (p*nvalues)/(ix+1)        
        p_adj_list.append(fdr_adj_p_val)
    
    return p_adj_list