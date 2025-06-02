import pandas as pd
from sli_algo import SliAlgo
from datetime import datetime

def main():
    # This is the list of KO genes that are to be tested for each mutant, it contains 10053 rows
    ko_genes_to_include = pd.read_excel('data/interim/reactome genes list (to include).xlsx')['gene_name'].to_list()

    # This is the list of mutant genes to pass to the algorithm
    # mutants_to_include = pd.read_csv("sli-algo inputs/unique_mutants_tested.csv", header=None)
    # mutants_to_include = mutants_to_include.iloc[:,0].to_list()
    mutants_to_include = ["ARID1A"]

    data_path = "data/interim/final_table_newdata.csv"

    algo = SliAlgo(data_path, mutants_to_include, ko_genes_to_include)

    print("Loading data...")
    algo.load_data()

    # Run the algorithm
    results_df = algo.run_crispr_database(low_percentile=10, high_percentile=90, threads = 6)

    # Save results to Excel file
    save_path = "data/processed/sli_algo_output"
    try:
        results_df.to_csv(save_path)

    except FileNotFoundError:
        print(f"Warning: unknown file path {save_path} \n Dumped output.csv at script location.")
        results_df.to_csv("output.csv")

    except Exception as e:
        print("Unable to save file: ", e)

if __name__=="__main__":
    main()