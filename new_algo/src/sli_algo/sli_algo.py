import pandas as pd
import numpy as np
import multiprocessing
from scipy.stats import mannwhitneyu 
from statsmodels.stats.multitest import multipletests
from typing import List

class SliAlgo:
    def __init__(
        self, 
        filtered_data_path:str, 
        mutants_to_include:List[str], 
        ko_genes_to_include:List[str]
    ):
        self.filtered_data_path = filtered_data_path
        self.mutants_to_include = mutants_to_include
        self.ko_genes_to_include = ko_genes_to_include

    def load_data(self):
        try:
            self.filtered_data = pd.read_csv(self.filtered_data_path)
        except FileNotFoundError:
            raise ValueError(f"File for Final Table not found: {self.filtered_data_path}")
        # Preprocess data into a dictionary for quick access
        self.gene_groups = {gene: group for gene, group in self.filtered_data.groupby('gene_name')}


    def split_populations(self, mutant_data, lower_percentile, upper_percentile):
        gene_expr = mutant_data['gene_expression'].values
        low_bound = np.percentile(gene_expr, lower_percentile)
        high_bound = np.percentile(gene_expr, upper_percentile)
        
        mask_low = gene_expr <= low_bound
        mask_high = gene_expr > high_bound
        
        population = np.where(mask_low, 'low', np.where(mask_high, 'high', None))
        mutant_data = mutant_data.assign(population=population)
        return mutant_data.dropna(subset=['population'])

    def one_sided_test(self, mutant, gene, data):
        data = data.drop_duplicates(subset=["DepMap_ID", "gene_name", "crispr_effect", "population"])
        if data.empty or data["population"].nunique() < 2:
            return None
        
        counts = data["population"].value_counts()
        high_count = counts.get("high", 0)
        low_count = counts.get("low", 0)
        if high_count == 0 or low_count == 0:
            return None
        
        try:
            stat, p_value = mannwhitneyu(
                data[data["population"] == "low"]["crispr_effect"],
                data[data["population"] == "high"]["crispr_effect"],
                alternative='less',
                method='asymptotic'
            )
        except ValueError:
            return None
        
        result = {
            "mutant": mutant,
            "gene": gene,
            "high": high_count,
            "low": low_count,
            "p_value": p_value,
            "statistic": stat,
        }

        return result

    def run_hypothesis_test_unique_percentiles(
            self, 
            mutant:str, 
            gene:str, 
            low_percentile:int, 
            high_percentile:int
        ):
        # Fetch CRISPR data
        crispr_data = self.gene_groups[gene].dropna(subset=['crispr_effect'])
        if crispr_data.empty:
            return self.create_empty_result(mutant, gene, low_percentile, high_percentile)
        
        # Fetch mutant data 
        mutant_data = self.gene_groups.get(mutant, pd.DataFrame())
        # Filter it to eliminate entries that aren't in CRISPR_data 
        mutant_data = mutant_data[mutant_data['DepMap_ID'].isin(crispr_data['DepMap_ID'])].dropna(subset=['gene_expression'])
        if mutant_data.empty:
            return self.create_empty_result(mutant, gene, low_percentile, high_percentile)
        
        # 3. Split populations
        expression_populations = split_populations(mutant_data, low_percentile, high_percentile)
        if expression_populations.empty:
            return self.create_empty_result(mutant, gene, low_percentile, high_percentile)
        
        # 4. Map population to CRISPR data
        population_map = expression_populations.set_index('DepMap_ID')['population']
        crispr_population_data = crispr_data[crispr_data['DepMap_ID'].isin(population_map.index)].copy()
        crispr_population_data['population'] = crispr_population_data['DepMap_ID'].map(population_map)
        crispr_population_data.dropna(subset=['population'], inplace=True)
        
        # 5. Compute mean/median differences
        mean_median = crispr_population_data.groupby('population')['crispr_effect'].agg(['mean', 'median']).diff().fillna(0)
        diff_mean = mean_median['mean'].iloc[-1]
        diff_median = mean_median['median'].iloc[-1]
        
        # 6. Perform test
        result = self.one_sided_test(mutant, crispr_population_data, gene)
        if result is None:
            return self.create_empty_result(mutant, gene, low_percentile, high_percentile)
        
        result['low_percentile'] = low_percentile
        result['high_percentile'] = high_percentile
        result['diff_mean'] = diff_mean
        result['diff_median'] = diff_median
        return result

    def create_empty_result(self, mutant, gene, low, high):
        result = {
            "mutant": [mutant],
            "gene": [gene],
            "high": [0],
            "low": [0],
            "p_value": [np.nan],
            "statistic": [np.nan],
            "low_percentile": [low],
            "high_percentile": [high],
            "diff_mean": [0],
            "diff_median": [0]
        }

        return result

    # Save results to Excel file with one tab for each mutant
    def save_results_to_excel(self, results_df, filename="sli-algo outputs/results.xlsx"):
        # Create Excel writer
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            # Create a sheet with all results
            # First convert adjusted_p_value to float for proper sorting
            if 'adjusted_p_value' in results_df.columns:
                # Create a temporary float column for sorting
                results_df['adjusted_p_value_float'] = results_df['adjusted_p_value'].apply(lambda x: float(x) if pd.notna(x) else float('inf'))
                # Sort the entire dataframe by mutant and then by p-value
                all_results = results_df.sort_values(by=['mutant', 'adjusted_p_value_float'])
                # Remove the float column
                all_results = all_results.drop(columns=['adjusted_p_value_float'])
            else:
                all_results = results_df
                
            # Save all results to first sheet
            all_results.to_excel(writer, sheet_name='All_Results', index=False)
            
            # Create separate sheets for each mutant
            for mutant in results_df['mutant'].unique():
                # Filter data for this mutant
                mutant_data = results_df[results_df['mutant'] == mutant]
                
                # Sort by p-value numerically (ascending)
                if 'adjusted_p_value' in mutant_data.columns:
                    # Create a temporary float column for sorting
                    mutant_data['adjusted_p_value_float'] = mutant_data['adjusted_p_value'].apply(
                        lambda x: float(x) if pd.notna(x) else float('inf'))
                    # Sort by the numerical p-value
                    mutant_data = mutant_data.sort_values(by='adjusted_p_value_float')
                    # Remove the temporary float column
                    mutant_data = mutant_data.drop(columns=['adjusted_p_value_float'])
                
                # Save to sheet named with the mutant's name
                mutant_data.to_excel(writer, sheet_name=f'{mutant}', index=False)
            
            # Fix openpyxl bug: ensure all sheets are marked as visible before saving
            for ws in writer.book.worksheets:
                ws.sheet_state = 'visible'
                
        print(f"Results saved to {filename}")

    # Define the process_row function at module level for multiprocessing to work
    def process_row(self, row):
        return self.run_hypothesis_test_unique_percentiles(
            row["mutant"], row["gene"], row["low_percentile"], row["high_percentile"]
        )

    def run_crispr_database(
            self, 
            low_percentile : int, 
            high_percentile : int, 
            threads: int = 1
        ):  
        # Create a list of all SL pairs to be tested (excluding self-pairs)
        mutant_gene_pairs = [(mutant, ko_gene, low_percentile, high_percentile)
                for mutant in self.mutants_to_include
                for ko_gene in self.ko_genes_to_include if ko_gene != mutant]
        
        # Convert the list of tuples into a DataFrame
        mutant_gene_pairs_df = pd.DataFrame(mutant_gene_pairs, columns=["mutant", "gene", "low_percentile", "high_percentile"])

        rows = mutant_gene_pairs_df.to_dict('records')

        # Single-threaded execution
        if threads == 1:

            results_list = map(self.process_row, rows)

            # Convert the list of results to a DataFrame
            results_df = pd.DataFrame(results_list)
        
        elif threads > 1:
            # Convert DataFrame to list of dictionaries for multiprocessing
        
            # Create a pool of workers
            with multiprocessing.Pool(processes=threads) as pool:
                # Map the function to the rows
                results_list = pool.map(self.process_row, rows)
            
            # Convert the list of results to a DataFramep
            results_df = pd.DataFrame(results_list)
        
        else:
            raise ValueError("Invalid number of threads")

        # Apply Benjamini-Hochberg correction to p-values
        for mutant in self.mutants_to_include: 
            p_values = results_df.loc[results_df['mutant'] == mutant, 'p_value'].values
            reject, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
            results_df.loc[results_df['mutant'] == mutant, 'adjusted_p_value'] = corrected_p_values

        # Format adjusted p-values in scientific notation
        results_df['adjusted_p_value'] = results_df['adjusted_p_value'].apply(lambda x: f"{x:.6E}")

        return results_df


