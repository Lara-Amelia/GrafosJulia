import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
import numpy as np

def prepare_comparison_data(files_dict):
    all_dfs = []
    for label, path in files_dict.items():
        try:
            df = pd.read_csv(path)
            df.columns = df.columns.str.strip()
            df['Method'] = label
            
            def extract_metadata(inst):
                # Handles 'bi_a100_b500' -> size 600
                size_match = re.search(r'a(\d+)_b(\d+)', inst)
                if size_match:
                    n = int(size_match.group(1)) + int(size_match.group(2))
                else:
                    n_match = re.search(r'n(\d+)', inst)
                    n = int(n_match.group(1)) if n_match else 0
                
                # Handles 'p01%' -> 1.0
                p_match = re.search(r'p(\d+)%', inst)
                p = int(p_match.group(1)) if p_match else 0
                return pd.Series([n, p])

            df[['Size', 'Prob']] = df['instancia'].apply(extract_metadata)
            all_dfs.append(df[['Method', 'Size', 'Prob', 'mean_chi', 'mean_time']])
        except FileNotFoundError:
            print(f"Note: {path} not found. Skipping {label} for now.")
            continue
    
    return pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()

# 1. Define your comparison methods
comparison_methods = {
    "Standard GA (60%)": "results_GA_Bi_pmutation06_stag20_progresso1.csv",
    "GABinomial (80%)": "results_GABinomial_Bi_pmutation08_stag20_progresso1.csv",
    "BRKGA (Mutation)": "results_BRKGAmutation_Bi_stag20_progresso1.csv",
    "Pure BRKGA": "results_BRKGA_Bi_stag20_progresso1.csv"
}

# 2. Load and process main data
df_plot = prepare_comparison_data(comparison_methods)

# 3. Plotting logic
if not df_plot.empty:
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    sns.set_theme(style="whitegrid")

    # --- PLOT 1: QUALITY VS SIZE (LINEAR ZOOMED) ---
    sns.lineplot(data=df_plot, x="Size", y="mean_chi", hue="Method", marker="o", ax=axes[0])
    
    # Zoom logic: Focus on the metaheuristic range
    y_max = df_plot['mean_chi'].max() * 1.1 
    axes[0].set_ylim(0, y_max) 
    axes[0].set_title("Coloring Quality vs Size")
    axes[0].set_ylabel(r"Mean Chi ($\chi$)")

    # --- PLOT 2: QUALITY VS DENSITY (LINEAR ZOOMED) ---
    # FIXED: Using df_plot consistently to avoid NameError
    sns.lineplot(data=df_plot, x="Prob", y="mean_chi", hue="Method", marker="s", ax=axes[1])
    
    axes[1].set_ylim(0, y_max)
    axes[1].set_title("Coloring Quality vs Density (%)")
    axes[1].set_ylabel(r"Mean Chi ($\chi$)")

    plt.tight_layout()
    plt.savefig("bipartite_comparison_final.png", dpi=300)
    print("Comparison plots generated as 'bipartite_comparison_final.png'.")
else:
    print("Error: No data loaded. Check your CSV filenames.")