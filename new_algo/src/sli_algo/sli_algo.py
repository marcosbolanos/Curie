import pandas as pd
import numpy as np
import multiprocessing
import functools
import json
import traceback
from scipy.stats import mannwhitneyu 
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from datetime import datetime
from typing import List

class SliProcessor:
    """Lightweight processor that contains only the computation logic"""
    
    def __init__(self, gene_groups: dict):
        self.gene_groups = gene_groups
    
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

    def create_empty_result(self, mutant, gene, low, high):
        result = {
            "mutant": mutant,
            "gene": gene,
            "high": 0,
            "low": 0,
            "p_value": np.nan,
            "statistic": np.nan,
            "low_percentile": low,
            "high_percentile": high,
            "diff_mean": 0,
            "diff_median": 0
        }
        return result

    def run_hypothesis_test_unique_percentiles(
            self, 
            mutant: str, 
            gene: str, 
            low_percentile: int, 
            high_percentile: int
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
        expression_populations = self.split_populations(mutant_data, low_percentile, high_percentile)
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
        result = self.one_sided_test(mutant, gene, crispr_population_data)
        if result is None:
            return self.create_empty_result(mutant, gene, low_percentile, high_percentile)
        
        result['low_percentile'] = low_percentile
        result['high_percentile'] = high_percentile
        result['diff_mean'] = diff_mean
        result['diff_median'] = diff_median
        return result


def _process_row_multiprocessing(filtered_data_path: str, row: dict):
    """Module-level function for multiprocessing that loads data efficiently"""
    try:
        # Load data once per process (this gets cached by the OS)
        filtered_data = pd.read_csv(filtered_data_path)
        gene_groups = {gene: group for gene, group in filtered_data.groupby('gene_name')}
        
        # Create processor and run computation
        processor = SliProcessor(gene_groups)
        return processor.run_hypothesis_test_unique_percentiles(
            row["mutant"], row["gene"], row["low_percentile"], row["high_percentile"]
        )
    except Exception as e:
        error_info = {
            "row": row,
            "error": str(e),
            "traceback": traceback.format_exc(),
            "timestamp": datetime.now().isoformat()
        }
        return {"error": error_info}

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
        # Create processor instance
        self.processor = SliProcessor(self.gene_groups)


    # Delegate computation methods to processo
    def split_populations(self, mutant_data, lower_percentile, upper_percentile):
        return self.processor.split_populations(mutant_data, lower_percentile, upper_percentile)

    def one_sided_test(self, mutant, gene, data):
        return self.processor.one_sided_test(mutant, gene, data)

    def run_hypothesis_test_unique_percentiles(self, mutant: str, gene: str, low_percentile: int, high_percentile: int):
        return self.processor.run_hypothesis_test_unique_percentiles(mutant, gene, low_percentile, high_percentile)

    def create_empty_result(self, mutant, gene, low, high):
        return self.processor.create_empty_result(mutant, gene, low, high)

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
        try:
            return run_hypothesis_test_unique_percentiles(
                row["mutant"], row["gene"], row["low_percentile"], row["high_percentile"]
            )
        except Exception as e:
            error_info = {
                "row": row,
                "error": str(e),
                "traceback": traceback.format_exc(),
                "timestamp": datetime.now().isoformat()
            }
            error_data = {"error": error_info}
            print(error_data)
            return (error_data)

    def run_crispr_database(
            self, 
            low_percentile: int, 
            high_percentile: int, 
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
            print(f"Processing {len(rows)} gene pairs with single thread...")
            results_list = []
            for row in tqdm(rows, desc="Processing"):
                try:
                    result = self.run_hypothesis_test_unique_percentiles(
                        row["mutant"], row["gene"], row["low_percentile"], row["high_percentile"]
                    )
                    results_list.append(result)
                except Exception as e:
                    error_info = {
                        "row": row,
                        "error": str(e),
                        "traceback": traceback.format_exc(),
                        "timestamp": datetime.now().isoformat()
                    }
                    results_list.append({"error": error_info})
        
        elif threads > 1:
            print(f"Processing {len(rows)} gene pairs with {threads} threads...")
            # Use partial function with only the data path for memory efficiency
            process_func = functools.partial(_process_row_multiprocessing, self.filtered_data_path)
            
            with multiprocessing.Pool(processes=threads) as pool:
                results_list = list(tqdm(pool.imap(process_func, rows), total=len(rows), desc="Processing"))

        else:
            raise ValueError("Invalid number of threads")
        
        # Separate successful results from errors
        successful_results = [r for r in results_list if "error" not in r]
        errors = [r["error"] for r in results_list if "error" in r]

        # Print summary statistics
        total_processed = len(results_list)
        successful_count = len(successful_results)
        error_count = len(errors)
        
        print(f"\nProcessing Summary:")
        print(f"  Total processed: {total_processed}")
        print(f"  Successful: {successful_count}")
        print(f"  Errors: {error_count}")
        print(f"  Success rate: {(successful_count/total_processed)*100:.2f}%")

        # Log errors if any occurred
        if errors:
            error_filename = f"{self.filtered_data_path}_errors.json"
            try:
                # Try to load existing errors
                try:
                    with open(error_filename, 'r') as f:
                        existing_errors = json.load(f)
                except FileNotFoundError:
                    existing_errors = []
                
                # Append new errors
                existing_errors.extend(errors)
                
                # Save all errors
                with open(error_filename, 'w') as f:
                    json.dump(existing_errors, f, indent=2)
                
                print(f"Warning: {error_count} errors occurred during execution. Details saved to {error_filename}")
                
            except Exception as json_error:
                print(f"Warning: Could not save errors to JSON file: {json_error}")
                print(f"Errors that occurred: {error_count}")
        else:
            print("No errors occurred during processing.")
        
        # Convert the list of successful results to a DataFrame
        results_df = pd.DataFrame(successful_results) if successful_results else pd.DataFrame()

        # Apply Benjamini-Hochberg correction to p-values
        for mutant in mutants_to_include: 
            p_values = results_df.loc[results_df['mutant'] == mutant, 'p_value'].values
            if len(p_values) == 0:
                results_df.loc[results_df['mutant'] == mutant, 'adjusted_p_value'] = 'NA'
                print(f"Warning : no p-values found for mutant: {mutant}")
                continue
            reject, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
            results_df.loc[results_df['mutant'] == mutant, 'adjusted_p_value'] = corrected_p_values

        # Format adjusted p-values in scientific notation
        results_df['adjusted_p_value'] = results_df['adjusted_p_value'].apply(lambda x: f"{x:.6E}")

        return results_df


