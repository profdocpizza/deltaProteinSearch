import os
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
from collections import Counter
from wordcloud import WordCloud
import seaborn as sns
import matplotlib.colors as mcolors


def plot_cath_sunburst(df, output_dir):
    """
    Generates a sunburst plot for the CATH classification hierarchy.

    It splits the 'cath' column into hierarchical levels, groups by those levels,
    and creates a Plotly sunburst plot.
    """
    # Preprocess CATH codes into hierarchical levels if not already done
    if "cath_levels" not in df.columns:
        df["cath_levels"] = df["cath"].str.split(".")
        df[["Class", "Architecture", "Topology", "Superfamily"]] = pd.DataFrame(
            df["cath_levels"].tolist(), index=df.index
        )

    # Prepare data for sunburst
    cath_hierarchy = (
        df.groupby(["Class", "Architecture", "Topology", "Superfamily"])
        .size()
        .reset_index(name="Count")
    )

    # Create the sunburst plot with Plotly Express
    fig = px.sunburst(
        cath_hierarchy,
        path=["Class", "Architecture", "Topology", "Superfamily"],
        values="Count",
        title="CATH Classification Hierarchy",
        color="Class",
        color_discrete_sequence=px.colors.qualitative.Set2,
        width=800,
        height=800,
    )
    fig.update_layout(title_font_size=16, margin=dict(t=50, l=25, r=25, b=25))
    # Save as PNG and SVG
    for ext in ["png", "svg"]:
        output_path = os.path.join(output_dir, f"cath_sunburst.{ext}")
        fig.write_image(output_path, format=ext, scale=2)
    # fig.show()


def plot_ec_sunburst(df, output_dir):
    """
    Generates a sunburst plot for EC numbers by splitting them into hierarchical levels.

    This function splits the full EC number (dot-delimited) into four levels,
    padding missing parts with "Missing", then groups the data and plots the sunburst.
    """

    def split_ec(ec):
        # If EC number is missing, return all levels as "Missing"
        if pd.isna(ec):
            return ["Missing", "Missing", "Missing", "Missing"]
        parts = str(ec).split(".")
        while len(parts) < 4:
            parts.append("Missing")
        return parts[:4]

    # Create a new column with split EC levels
    df["EC_levels"] = df["ec_number"].apply(split_ec)
    df[["EC_Level1", "EC_Level2", "EC_Level3", "EC_Level4"]] = pd.DataFrame(
        df["EC_levels"].tolist(), index=df.index
    )

    # Prepare data for sunburst
    ec_hierarchy = (
        df.groupby(["EC_Level1", "EC_Level2", "EC_Level3", "EC_Level4"])
        .size()
        .reset_index(name="Count")
    )

    # Create the sunburst plot with Plotly Express
    fig = px.sunburst(
        ec_hierarchy,
        path=["EC_Level1", "EC_Level2", "EC_Level3", "EC_Level4"],
        values="Count",
        title="EC Classification Hierarchy",
        color="EC_Level1",
        color_discrete_sequence=px.colors.qualitative.Vivid,
        width=800,
        height=800,
    )
    fig.update_layout(title_font_size=16, margin=dict(t=50, l=25, r=25, b=25))

    # Save as PNG and SVG
    for ext in ["png", "svg"]:
        output_path = os.path.join(output_dir, f"ec_sunburst.{ext}")
        fig.write_image(output_path, format=ext, scale=2)

    # fig.show()


# def plot_keyword_wordcloud(df, output_dir):
#     """
#     Generates a keyword word cloud from the 'keywords' column and saves it as PNG and SVG.

#     The function:
#       - Splits the keywords string into individual keywords.
#       - Uses a Counter to compute keyword frequencies.
#       - Generates a word cloud using the WordCloud library.
#       - Saves the plot in both PNG and SVG formats in the specified output directory.
#     """
#     # Ensure the output directory exists
#     os.makedirs(output_dir, exist_ok=True)

#     # Process keywords: drop missing values, split by '; ', explode the list, and strip whitespace
#     keywords_series = df["keywords"].dropna().str.split("; ").explode().str.strip()
#     # remove keywords in a list
#     keywords_series = keywords_series[
#         ~keywords_series.isin(
#             [
#                 "3D-structure",
#                 "Reference proteome",
#                 "Direct protein sequencing",
#                 "Proteomics identification",
#             ]
#         )
#     ]
#     keyword_counts = Counter(keywords_series)

#     # Generate word cloud from frequencies
#     wc = WordCloud(
#         width=800,
#         height=400,
#         background_color="white",
#         colormap="plasma",
#         min_font_size=10,
#     )
#     wc_image = wc.generate_from_frequencies(keyword_counts)

#     plt.figure(figsize=(10, 6))
#     plt.imshow(wc_image, interpolation="bilinear")
#     plt.axis("off")
#     # plt.title("Keyword Word Cloud", fontsize=14, pad=10)
#     plt.tight_layout()


#     # Save as PNG and SVG
#     for ext in ["png", "svg"]:
#         output_path = os.path.join(output_dir, f"keyword_wordcloud.{ext}")
#         plt.savefig(output_path, format=ext, dpi=300)
#     # plt.show()
def plot_keyword_wordcloud(df, output_dir):
    """
    Generates a keyword word cloud from the 'keywords' column and saves it as PNG and SVG.

    The function:
      - Reads the keyword mapping (kw_df) to map each keyword to its broader category.
      - Splits the keywords string into individual keywords.
      - Removes a few specified keywords.
      - Uses a custom color function to assign each word a color based on its category
        using a colorblind-friendly palette.
      - Generates the word cloud using the WordCloud library.
      - Saves the plot in both PNG and SVG formats in the specified output directory.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read the keyword mapping file and rename columns as needed
    kw_df = pd.read_csv(
        "/home/tadas/code/deltaProteinSearch/uniprot_keywords.tsv", sep="\t"
    )
    kw_df.rename(
        columns={"Name": "keyword", "Category": "keyword_category"}, inplace=True
    )

    # Check if 'keywords' column exists
    if "keywords" not in df.columns:
        print("The DataFrame does not have a 'keywords' column.")
        return

    # Process keywords: drop missing values, split by '; ', explode, and strip whitespace
    keywords_series = df["keywords"].dropna().str.split("; ").explode().str.strip()

    # Remove some unwanted keywords
    remove_keywords = [
        "3D-structure",
        "Reference proteome",
        "Direct protein sequencing",
        "Proteomics identification",
    ]
    keywords_series = keywords_series[~keywords_series.isin(remove_keywords)]

    # Check if any keywords remain after filtering
    if keywords_series.empty:
        print("No keywords available after filtering. Word cloud cannot be generated.")
        return

    # Create a DataFrame from the keywords
    keywords_df = pd.DataFrame(keywords_series, columns=["keyword"])

    # Merge with the keyword mapping DataFrame to get the broader category for each keyword
    merged = keywords_df.merge(kw_df, how="left", on="keyword")
    # Fill missing keyword_category values with "Unknown"
    merged["keyword_category"] = merged["keyword_category"].fillna("Unknown")

    # Count the frequency of each keyword
    keyword_counts = merged["keyword"].value_counts().to_dict()

    # If no keywords are found, warn and exit
    if not keyword_counts:
        print("No keywords found after processing. Word cloud cannot be generated.")
        return

    # Create a mapping from keyword to its category (each keyword should belong to the same category)
    keyword_to_category = (
        merged.drop_duplicates("keyword")
        .set_index("keyword")["keyword_category"]
        .to_dict()
    )

    # Get the unique categories and assign each a color from a colorblind-friendly palette
    unique_categories = sorted(merged["keyword_category"].unique())
    num_categories = len(unique_categories)
    palette = sns.color_palette("colorblind", n_colors=num_categories)
    # Convert palette colors to hexadecimal strings
    category_color_map = {
        cat: mcolors.to_hex(color) for cat, color in zip(unique_categories, palette)
    }

    # Define a custom color function for WordCloud
    def custom_color_func(
        word, font_size, position, orientation, random_state=None, **kwargs
    ):
        # Look up the keyword's category; default to "Unknown" if not found
        cat = keyword_to_category.get(word, "Unknown")
        # Return the corresponding color (defaulting to black if category not mapped)
        return category_color_map.get(cat, "#000000")

    # Generate the word cloud using the custom color function
    wc = WordCloud(
        width=800,
        height=400,
        background_color="white",
        min_font_size=10,
        color_func=custom_color_func,
    )
    wc_image = wc.generate_from_frequencies(keyword_counts)

    # Plot the word cloud
    plt.figure(figsize=(10, 6))
    plt.imshow(wc_image, interpolation="bilinear")
    plt.axis("off")
    plt.tight_layout()

    # Save the word cloud as PNG and SVG
    for ext in ["png", "svg"]:
        output_path = os.path.join(output_dir, f"keyword_wordcloud.{ext}")
        plt.savefig(output_path, format=ext, dpi=300)


def plot_keyword_stacked_bar(df, output_dir):
    kw_df = pd.read_csv(
        "/home/tadas/code/deltaProteinSearch/uniprot_keywords.tsv", sep="\t"
    )
    kw_df.rename(
        columns={"Name": "keyword", "Category": "keyword_category"}, inplace=True
    )

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process the keywords: drop rows with missing keywords, split by '; ', explode, and strip whitespace
    df_keywords = df[["keywords"]].dropna().copy()
    df_keywords["keyword"] = df_keywords["keywords"].str.split("; ")
    df_keywords = df_keywords.explode("keyword")
    df_keywords["keyword"] = df_keywords["keyword"].str.strip()

    # Merge with the keyword mapping DataFrame (kw_df)
    merged_df = df_keywords.merge(kw_df, how="left", on="keyword")

    # Fill missing keyword_category values with "Unknown"
    merged_df["keyword_category"] = merged_df["keyword_category"].fillna("Unknown")

    # Group by 'keyword_category' and 'keyword' to count occurrences
    keyword_group = (
        merged_df.groupby(["keyword_category", "keyword"])
        .size()
        .reset_index(name="Count")
    )

    # Pivot the data so that each row corresponds to a keyword_category
    # and each column represents an individual keyword
    keyword_pivot = keyword_group.pivot(
        index="keyword_category", columns="keyword", values="Count"
    ).fillna(0)

    # Plot the stacked bar chart using Matplotlib (pandas built-in plotting)
    ax = keyword_pivot.plot(
        kind="bar", stacked=True, figsize=(12, 8), colormap="viridis"
    )
    plt.title("Stacked Bar Chart of Keywords by Broader Category", fontsize=14, pad=10)
    plt.xlabel("Keyword Category", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # Save the chart as PNG and SVG
    for ext in ["png", "svg"]:
        output_path = os.path.join(output_dir, f"keyword_stacked_bar.{ext}")
        plt.savefig(output_path, format=ext, dpi=300)

    plt.show()


# Example usage
if __name__ == "__main__":
    # Read the dataset
    df = pd.read_csv(
        "/home/tadas/code/deltaProteinSearch/outputs/pdb_with_uniprot_details_async.csv"
    )

    # Specify the output directory for saving images
    output_directory = "/home/tadas/code/deltaProteinSearch/outputs/figures"
    # plot_keyword_stacked_bar(df, output_directory)
    # Plot and show CATH sunburst visualization
    plot_cath_sunburst(df, output_directory)

    # Plot and show EC sunburst visualization
    plot_ec_sunburst(df, output_directory)

    # Generate, display, and save the keyword word cloud
    plot_keyword_wordcloud(df, output_directory)
