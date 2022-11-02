import pandas as pd
import plotly.express as px
from bioinfokit.analys import stat
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statistical.fdr import get_false_discovery_rate

def anova(df: pd.DataFrame, column_classes: str) -> pd.DataFrame:
    """
    _summary_

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    column_classes : str
        _description_

    Returns
    -------
    pd.DataFrame
        _description_
    """
    
    res = stat()
    
    anovastats = {}
    for x in df:
        model = "{} ~ {}".format(x, column_classes)
        if x != column_classes:
            res.anova_stat(df=df, res_var= x , anova_model=model)
            anovastats[x] = res.anova_summary["PR(>F)"][0]

    print(anovastats)

    pval = pd.DataFrame.from_dict(anovastats, orient='index')

    # Calculate FDR for each p-value
    pval["FDR"] = get_false_discovery_rate(pval, 0)
    
    return pval

def render_and_tukey_top_fdr(df: pd.DataFrame, column_classes: str, pval: pd.DataFrame, fdr_limit: float = 0.05):
    """
    _summary_

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    column_classes : str
        _description_
    pval : pd.DataFrame
        _description_
    fdr_limit : float, optional
        _description_, by default 0.05
    """
    padj005 = pval[pval['FDR'] <= fdr_limit]
    padj005.head()
    for names in padj005.index.values:
        tukey = pairwise_tukeyhsd(endog=df[names], groups=df[column_classes], alpha=0.05)
        print(names)
        print(tukey)

    for mass in padj005.index.values:
        fig = px.box(df, x= column_classes, y= mass, color= column_classes)
        fig.show()
