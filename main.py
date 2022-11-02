from load_data import load_dataframe_from_csv, is_binary_classification
from statistical.anova import anova, render_and_tukey_top_fdr
from statistical.ttest import ttest, render_top_fdr

def main():
    ## Asking file name
    file = input("Please, enter file name: ")    
    
    # Prepare dataframe
    df = load_dataframe_from_csv(file)

    ## Asking column name with variables
    column_classes = input("Please, enter column name with variables: ")
    if not is_binary_classification(df, column_classes):
        print("Performing Anova and Tukey tests")
        pval = anova(df, column_classes)
        render_and_tukey_top_fdr(df, column_classes, pval)

    else:
        print("Performing T-test")
        pval = ttest(df, column_classes)
        render_top_fdr(df, column_classes, pval)

if __name__ == "__main__":
    main()