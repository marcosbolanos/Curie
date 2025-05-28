import pandas as pd
from src.sli_algo import SliAlgo

def main():
    # this is the list of KO genes that are to be tested for each mutant, it contains 10053 rows
    ko_genes_to_include = pd.read_excel('data/interim/reactome genes list (to include).xlsx')['gene_name'].to_list()

    # mutants_to_include = pd.read_csv("sli-algo inputs/unique_mutants_tested.csv", header=None)
    # mutants_to_include = mutants_to_include.iloc[:,0].to_list()

    mutants_to_include = ["ARID1A"]

    data_path = "data/interim/final_table_newdata.csv"

    algo = SliAlgo(data_path, mutants_to_include, ko_genes_to_include)
    algo.load_data()
    results_df = algo.run_hypothesis_test_unique_percentiles("SMARCA2", "SMARCA4", 10, 90)
    print(results_df)
    # results_df = algo.run_crispr_database(low_percentile=10, high_percentile=90, threads = 6)

    # Save results to Excel file
    filename = "test.csv"
    results_df.to_csv(filename)

if __name__=="__main__":
    main()