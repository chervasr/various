import pandas as pd

def load_dataframe_from_csv(file: str):
    df: pd.DataFrame = pd.read_csv(file)
    df = df.transpose()
    df = df.rename(columns=df.iloc[0])
    df = df.drop(["Mass"],axis=0)
    df = df.astype({col: float for col in df.columns[1:]})
    return df

def is_binary_classification(df, column_name):
    return len(df[column_name].value_counts(ascending=True)) == 2