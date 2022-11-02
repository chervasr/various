import pandas as pd
import plotly.express as px
import scipy
from statistical.fdr import get_false_discovery_rate

def ttest(df: pd.DataFrame, column_classes: str) -> pd.DataFrame:
    """
    Perform T-test for two groups after a binary classification,for all columns in a dataframe of mass spectometry. It also calculates False Discovery Rate for the pvalues of this columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a peak matrix from mass spectometry. Each column represent a different peak of the graph, with one column containing class names. Each row corresponds with a sample.
    column_classes : str
        Column name of the column that difference the samples into two classes.

    Returns
    -------
    pd.DataFrame
        DataFrame containing P-Value and FDR for each peak.
    """
    # Binary mapper for classes
    ix_a = df[column_classes] == df[column_classes][1]
    
    # Calculate p-value for each mass
    tstats = {}
    for x in df:
        if x != column_classes:
            tstats[x] = scipy.stats.ttest_ind(df[x][ix_a], df[x][~ix_a]).pvalue
    pval = pd.DataFrame.from_dict(tstats, orient='index')
    pval.rename(columns={ pval.columns[0]: "P-Value" }, inplace = True)

    # Calculate FDR for each p-value
    pval["FDR"] = get_false_discovery_rate(pval, "P-Value")
    
    return pval

def render_top_fdr(df: pd.DataFrame, column_classes: str, pval: pd.DataFrame, fdr_limit: float = 0.05):
    """
    Represents a boxplot for each of the peaks with a FDR lower than the value introducted, for both classes.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a peak matrix from mass spectometry. Each column represent a different peak of the graph, with one column containing class names. Each row corresponds with a sample.
    column_classes : str
        Column name of the column that difference the samples into two classes.
    pval : pd.DataFrame
        DataFrame with p-values and FDR for each peak
    fdr_limit : float, optional
        Maximun value of FDR to take in consideration, by default 0.05
    """
    padj005 = pval[pval['FDR'] <= fdr_limit]
    print(padj005)
    for mass in padj005.index.values:
        fig = px.box(df, x= column_classes, y= mass, color= column_classes)
        fig.show()
