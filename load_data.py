import pandas as pd

def load_dataframe_from_csv(file: str):
    """
    Load the content of a CSV file into a Pandas Dataframe. Adjust the Dataframe to the expected input

    Parameters
    ----------
    file : str
        Path to a CSV file

    Returns
    -------
    pd.DataFrame
        A Dataframe with where each row is ... and each column is ...
    """
    df: pd.DataFrame = pd.read_csv(file)
    df = df.transpose()
    df = df.rename(columns=df.iloc[0])
    df = df.drop(["Mass"],axis=0)
    df = df.astype({col: float for col in df.columns[1:]})
    return df

def is_binary_classification(df, column_name):
    """
    _summary_

    Parameters
    ----------
    df : _type_
        _description_
    column_name : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    return len(df[column_name].value_counts(ascending=True)) == 2
    # Alternative: df[column_name].nunique() == 2    
    # TODO What if column_name does not exist in df?
    # TODO: What if only one class?