""" LIBRARIES """
import pandas as pd
import numpy as np
import re
import glob
from collections import Counter


""" FUNCTIONS """
def load_data(file_path):
    columns_needed = ["chr", "pos", "ref", "alt", "rs_dbSNP151", "ANNOVAR_refseq_Gene_ID",
                      "ANNOVAR_refseq_Closest_gene(intergenic_only)", "SnpEff_refseq_Gene_ID", "VEP_refseq_Gene_Name"]
    df = pd.read_csv(file_path, usecols=columns_needed, delimiter="\t", dtype=object, compression="gzip")
    return df

def extract_refseq(lists):
    """
    Process RefSeq genes from input lists, filter, and polish them by extracting gene identifiers.

    Parameters
    ----------
    lists : list
        List of strings containing RefSeq gene information.

    Returns
    -------
    list
        List of processed RefSeq genes for each SNP.
    """
    processed_genes = []

    for entry in lists:
        # Handle empty or NaN entries
        if not entry or entry.lower() == "nan":
            processed_genes.append([])
            continue

        # Remove "CHR_START-" ONLY if it's at the beginning of the string
        entry = re.sub(r"^CHR_START-?", "", entry)

        # Split entry by both "|" and "," to handle multiple gene listings
        sub_entries = re.split(r"[|,]", entry)

        final_genes = []
        for sub_entry in sub_entries:
            # Remove parts that contain "NONE" to exclude non-informative entries
            if "NONE" in sub_entry.upper():
                continue

            # Extract only the gene name before ':' if present
            if ':' in sub_entry:
                sub_entry = sub_entry.split(':')[0]

            # Split by '-' but preserve valid numeric parts (e.g., "NKX2-3", "MIR3648-1")
            parts = sub_entry.split('-')
            final_parts = []
            current_part = parts[0]

            for i in range(1, len(parts)):
                # If the current part ends with a digit and the next part starts with a digit, merge them
                if current_part[-1].isdigit() and parts[i][0].isdigit():
                    current_part += "-" + parts[i]
                # If the next part is "AS" followed by a number, merge it back
                elif re.match(r"AS\d+", parts[i]):
                    current_part += "-" + parts[i]
                else:
                    final_parts.append(current_part)
                    current_part = parts[i]

            final_parts.append(current_part)  # Append the last part
            final_genes.extend(final_parts)

        # Remove empty strings and add to the final processed list
        processed_genes.append([gene for gene in final_genes if gene])

    return processed_genes

def load_gene_mapping(file_path):
    """
    Load gene mapping data from a file and create dictionaries to map:
    - Entrez ID → HGNC ID
    - Symbol → HGNC ID

    Parameters
    ----------
    file_path : str
        Path to the CSV file containing gene mapping data.

    Returns
    -------
    tuple of dicts
        (entrez_id → hgnc_id, symbol → hgnc_id)
    """
    # Load file with all columns as strings to prevent dtype issues
    df = pd.read_csv(file_path, dtype=str).fillna("")

    # Build Entrez-to-HGNC mapping (used for LOC-style entries)
    loc_mapping = {
        row["entrez_id"].strip(): row["uniprot_ids"].strip()
        for _, row in df.iterrows() if row["entrez_id"].strip()
    }

    # Build Symbol-to-HGNC mapping (used for standard gene names)
    symbol_mapping = {
        row["symbol"].strip(): row["uniprot_ids"].strip()
        for _, row in df.iterrows() if row["symbol"].strip()
    }

    return loc_mapping, symbol_mapping


def map_refseq_genes(refseq_genes_list, loc_mapping, symbol_mapping):
    """
    Map RefSeq gene names to Ensembl Gene IDs using the given mappings.

    Parameters
    ----------
    refseq_genes_list : list of lists
        Each sublist contains RefSeq gene names for a given SNP.
    loc_mapping : dict
        Mapping of LOC gene Entrez IDs to Ensembl Gene IDs.
    symbol_mapping : dict
        Mapping of gene symbols to Ensembl Gene IDs.

    Returns
    -------
    tuple:
        - list of lists: Mapped Ensembl Gene IDs per SNP (only successful matches).
        - int: Total number of unmatched gene IDs across all SNPs (SNP-level mismatches).
        - int: Number of unique unmatched gene IDs across the chromosome.
    """
    mapped_genes_list = []          # Final list to hold matched genes per SNP
    total_unmatched_count = 0       # Counter for unmatched gene IDs across all SNPs
    unmatched_genes_set = set()     # Set to track unique unmatched gene IDs

    # Process each SNP's gene list
    for gene_list in refseq_genes_list:
        matched = []                # Store matched gene IDs for this SNP
        seen_unmatched = set()      # Track unmatched gene IDs already counted for this SNP

        for gene in gene_list:
            mapped_gene = None

            # Check for LOC-prefixed gene and try to map using the Entrez ID
            if isinstance(gene, str) and gene.startswith("LOC") and gene[3:].isdigit():
                entrez_id = gene[3:]
                mapped_gene = loc_mapping.get(entrez_id)

            # If not a LOC match, try symbol mapping
            if not mapped_gene:
                mapped_gene = symbol_mapping.get(gene)

            # If mapping successful, keep it
            if mapped_gene:
                matched.append(mapped_gene)
            else:
                # Track unmatched genes for total SNP count (avoid counting duplicates per SNP)
                if gene not in seen_unmatched:
                    seen_unmatched.add(gene)
                    total_unmatched_count += 1

                # Track all unmatched genes for global uniqueness
                unmatched_genes_set.add(gene)

        mapped_genes_list.append(matched)

    # Return mapped results and mismatch tracking info
    return mapped_genes_list, total_unmatched_count, len(unmatched_genes_set)


def agreement_with_tools(AN, SN, VP, method="concatenated", track_over_annotation=False):
    """
    Track tool-specific agreement and ensure SNPs are uniquely categorized.

    Parameters
    ----------
    AN, SN, VP : list of lists
        Gene annotations from Annovar, SnpEff, and VEP respectively.
    method : str
        Agreement method ("concatenated", "majority_wins", or "combined").
    track_over_annotation : bool
        Whether to track tools that annotate beyond the reference set.

    Returns
    -------
    dict
        {
            agreement counts,
            over-annotation counts (if enabled),
            number of SNPs with no annotations at all
        }
    """
    # Initialize counters
    counters = {
        'A': 0,   # Annovar agrees alone
        'S': 0,   # SnpEff agrees alone
        'V': 0,   # VEP agrees alone
        'AS': 0,  # Annovar & SnpEff agree
        'AV': 0,  # Annovar & VEP agree
        'SV': 0,  # SnpEff & VEP agree
        'ASV': 0, # All tools agree
        'none': 0 # No tools agree
    }

    over_annotation = {'A': 0, 'S': 0, 'V': 0} if track_over_annotation else None
    empty_reference_count = 0

    # Ensure inputs are lists of lists
    AN = [list(item) if isinstance(item, (list, np.ndarray)) else [] for item in AN]
    SN = [list(item) if isinstance(item, (list, np.ndarray)) else [] for item in SN]
    VP = [list(item) if isinstance(item, (list, np.ndarray)) else [] for item in VP]

    # Iterate through SNPs
    for an, sn, vp in zip(AN, SN, VP):
        an, sn, vp = set(an), set(sn), set(vp)

        # Determine reference set
        if method == "concatenated":
            reference_set = an | sn | vp
        elif method in ["majority_wins", "combined"]:
            gene_counts = {}
            for gene in an | sn | vp:
                gene_counts[gene] = (gene in an) + (gene in sn) + (gene in vp)
            reference_set = {gene for gene, count in gene_counts.items() if count >= 2}
        else:
            raise ValueError("Invalid method. Use 'concatenated', 'majority_wins', or 'combined'.")

        if not reference_set:
            empty_reference_count += 1
            continue

        tools_agree = {
            'A': an == reference_set,
            'S': sn == reference_set,
            'V': vp == reference_set
        }

        if tools_agree['A'] and tools_agree['S'] and tools_agree['V']:
            counters['ASV'] += 1
        elif tools_agree['A'] and tools_agree['S']:
            counters['AS'] += 1
        elif tools_agree['A'] and tools_agree['V']:
            counters['AV'] += 1
        elif tools_agree['S'] and tools_agree['V']:
            counters['SV'] += 1
        elif tools_agree['A']:
            counters['A'] += 1
        elif tools_agree['S']:
            counters['S'] += 1
        elif tools_agree['V']:
            counters['V'] += 1
        else:
            counters['none'] += 1

        if track_over_annotation:
            for tool, annotations in zip(['A', 'S', 'V'], [an, sn, vp]):
                if not annotations.issubset(reference_set):
                    over_annotation[tool] += 1

    result = counters
    if track_over_annotation:
        result["Over-Annotation Counts"] = over_annotation
    result["Empty Reference Count"] = empty_reference_count

    return result



def agreement_rate(AN, SN, VP, method="concatenated", track_over_annotation=False):
    """
    Calculate the agreement counts of gene annotations per SNP for each tool,
    relative to the reference list, and optionally track over-annotation.
    Also tracks how many SNPs had an empty reference set.
    """
    counters = {
        'total_reference': 0,
        'agreement_counts': {'A': 0, 'S': 0, 'V': 0},
        'Empty Reference Count': 0
    }

    over_annotation = {'A': 0, 'S': 0, 'V': 0} if track_over_annotation else None

    # Clean input
    AN = [list(item) if isinstance(item, (list, np.ndarray)) else [] for item in AN]
    SN = [list(item) if isinstance(item, (list, np.ndarray)) else [] for item in SN]
    VP = [list(item) if isinstance(item, (list, np.ndarray)) else [] for item in VP]

    for an, sn, vp in zip(AN, SN, VP):
        an, sn, vp = set(an), set(sn), set(vp)

        # Determine the reference
        if method == "concatenated":
            reference_set = an | sn | vp
        elif method == "majority_wins":
            gene_counts = {}
            for gene in an | sn | vp:
                gene_counts[gene] = (gene in an) + (gene in sn) + (gene in vp)
            reference_set = {gene for gene, count in gene_counts.items() if count >= 2}
        else:
            raise ValueError("Invalid method. Use 'concatenated' or 'majority_wins'.")

        # Track empty reference
        if not reference_set:
            counters['Empty Reference Count'] += 1
            continue

        counters['total_reference'] += len(reference_set)

        # Check how many genes from reference each tool got
        for tool, annotations in zip(['A', 'S', 'V'], [an, sn, vp]):
            match_count = sum(1 for gene in reference_set if gene in annotations)
            counters['agreement_counts'][tool] += match_count

            if track_over_annotation:
                over_annotation[tool] += len(annotations - reference_set)

    result = {
        'Total Reference Size': counters['total_reference'],
        'Agreement Counts': counters['agreement_counts'],
        'Empty Reference Count': counters['Empty Reference Count']
    }

    if track_over_annotation:
        result['Over-Annotation Counts'] = over_annotation

    return result


def agreement_with_tools_sunburst(AN, SN, VP):
    """
    Tracks tool agreement for each SNP and generates labels for sunburst visualization.
    Uses the concatenated method (union of all tools) as the reference set.
    Returns a dictionary of counts.
    """
    counts = Counter()
    AN = [set(a) if isinstance(a, (list, np.ndarray)) else set() for a in AN]
    SN = [set(s) if isinstance(s, (list, np.ndarray)) else set() for s in SN]
    VP = [set(v) if isinstance(v, (list, np.ndarray)) else set() for v in VP]

    for an, sn, vp in zip(AN, SN, VP):
        ref_set = an | sn | vp
        if not ref_set:
            counts["Empty Reference Count"] += 1
            continue

        tools_agree = {
            'A': an == ref_set,
            'S': sn == ref_set,
            'V': vp == ref_set
        }

        subsets = {
            'A': bool(an) and an.issubset(ref_set),
            'S': bool(sn) and sn.issubset(ref_set),
            'V': bool(vp) and vp.issubset(ref_set)
        }

        full_agree = [k for k, v in tools_agree.items() if v]

        if len(full_agree) == 3:
            counts['ASV'] += 1
        elif len(full_agree) == 2:
            remaining = [k for k in ['A', 'S', 'V'] if k not in full_agree][0]
            label = "".join(sorted(full_agree)) + f"_with_{remaining}_subset" if subsets[remaining] else "".join(sorted(full_agree))
            counts[label] += 1
        elif len(full_agree) == 1:
            agreeing_tool = full_agree[0]
            others = [k for k in ['A', 'S', 'V'] if k != agreeing_tool and subsets[k]]
            if others:
                label = f"{agreeing_tool}_with_" + "_".join(sorted(others)) + "_subset"
                counts[label] += 1
            else:
                counts[agreeing_tool] += 1
        else:
            # No full match, but there might be partial overlaps
            overlap = (an & sn) | (an & vp) | (sn & vp)
            if overlap:
                counts["Partial_Agree_Overlap"] += 1
            else:
                counts["disjointed"] += 1

    return dict(counts)


def data_process(file_path, loc_mapping, symbol_mapping):
    """
    Process a single file and return results for genic and intergenic regions.
    """
    # Load data
    #file_path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\Code\\Test_data\\Test_data.txt.gz"
    #file_path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\AnnoQ_data\\21.annotated.snp.gz"
    cdata = load_data(file_path)

    # Load gene mappings (LOC-based and symbol-based) -  for debugging
    # loc_mapping, symbol_mapping = load_gene_mapping("C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\hgnc_complete_set.csv")

    # Define genic and intergenic regions based on proper classification
    AN_ID_genic = cdata["ANNOVAR_refseq_Gene_ID"]  # Gene IDs for genic regions
    AN_ID_intergenic = cdata["ANNOVAR_refseq_Closest_gene(intergenic_only)"]  # Closest gene for intergenic regions
    AN_ID_tog = [AN_ID_intergenic[i] if gene_id.startswith('.') else gene_id for i, gene_id in enumerate(AN_ID_genic)]

    SN_ID_tog = cdata["SnpEff_refseq_Gene_ID"]  # Gene IDs for SnpEff
    SN_ID_intergenic = [SN_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if gene_id.startswith('.')]
    SN_ID_genic = [SN_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]

    VP_ID_tog = cdata["VEP_refseq_Gene_Name"]  # Gene IDs for VEP
    VP_ID_intergenic = [VP_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if gene_id.startswith('.')]
    VP_ID_genic = [VP_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]

    # Fix genic and intergenic lists for ANNOVAR
    AN_ID_genic = [gene_id for gene_id in AN_ID_genic if not gene_id.startswith('.')]
    AN_ID_intergenic = [gene_id for gene_id in AN_ID_intergenic if not gene_id.startswith('.')]

    # ========== Extract and Map Gene IDs ==========
    AN_ID_genic, an_genic_unmatch, an_genic_unique = map_refseq_genes(extract_refseq(AN_ID_genic), loc_mapping, symbol_mapping)
    AN_ID_intergenic, an_inter_unmatch, an_inter_unique = map_refseq_genes(extract_refseq(AN_ID_intergenic), loc_mapping, symbol_mapping)

    SN_ID_genic, sn_genic_unmatch, sn_genic_unique = map_refseq_genes(extract_refseq(SN_ID_genic), loc_mapping, symbol_mapping)
    SN_ID_intergenic, sn_inter_unmatch, sn_inter_unique = map_refseq_genes(extract_refseq(SN_ID_intergenic), loc_mapping, symbol_mapping)

    VP_ID_genic, vp_genic_unmatch, vp_genic_unique = map_refseq_genes(extract_refseq(VP_ID_genic), loc_mapping, symbol_mapping)
    VP_ID_intergenic, vp_inter_unmatch, vp_inter_unique = map_refseq_genes(extract_refseq(VP_ID_intergenic), loc_mapping, symbol_mapping)


    # ========== Agreement Analyses (Concatenated Only) ==========
    genic_concat = agreement_with_tools(AN_ID_genic, SN_ID_genic, VP_ID_genic, method="concatenated")
    genic_concat_rate = agreement_rate(AN_ID_genic, SN_ID_genic, VP_ID_genic)

    inter_concat = agreement_with_tools(AN_ID_intergenic, SN_ID_intergenic, VP_ID_intergenic, method="concatenated")
    inter_concat_rate = agreement_rate(AN_ID_intergenic, SN_ID_intergenic, VP_ID_intergenic)

    # Sunburst breakdowns
    genic_sunburst = agreement_with_tools_sunburst(AN_ID_genic, SN_ID_genic, VP_ID_genic)
    intergenic_sunburst = agreement_with_tools_sunburst(AN_ID_intergenic, SN_ID_intergenic, VP_ID_intergenic)


    # Aggregate results
    chromosome = cdata["chr"].iloc[0]  # Assume all rows belong to the same chromosome

    return {
        "chromosome": chromosome,
        "genic_concat": genic_concat,
        "genic_concat_rate": genic_concat_rate,
        "genic_sunburst": genic_sunburst,
        "inter_concat": inter_concat,
        "inter_concat_rate": inter_concat_rate,
        "inter_sunburst": intergenic_sunburst,
        "genic_unmatched": {
            "an_genic_unmatch": an_genic_unmatch,
            "an_genic_unique": an_genic_unique,
            "sn_genic_unmatch": sn_genic_unmatch,
            "sn_genic_unique": sn_genic_unique,
            "vp_genic_unmatch": vp_genic_unmatch,
            "vp_genic_unique": vp_genic_unique,
        },
        "inter_unmatched": {
            "an_inter_unmatch": an_inter_unmatch,
            "an_inter_unique": an_inter_unique,
            "sn_inter_unmatch": sn_inter_unmatch,
            "sn_inter_unique": sn_inter_unique,
            "vp_inter_unmatch": vp_inter_unmatch,
            "vp_inter_unique": vp_inter_unique,
        }
    }


def flatten_results(data, result_type, region):
    """
    Flatten results for a specific analysis type and region into a DataFrame.

    Parameters
    ----------
    data : list
        List of dictionaries containing results for each chromosome.
    result_type : str
        Type of analysis (e.g., "concat", "concat_rate", "sunburst").
    region : str
        Region type ("genic" or "intergenic").

    Returns
    -------
    pd.DataFrame
        Flattened results as a DataFrame.
    """
    flattened = []
    for result in data:
        row = {"chromosome": result["chromosome"], "region": region}

        if result_type == "concat_rate":
            row.update(result[result_type]["Agreement Counts"])
            row["Total Reference Size"] = result[result_type]["Total Reference Size"]
            row["Empty Reference Count"] = result[result_type]["Empty Reference Count"]

        elif result_type == "sunburst":
            row.update(result[result_type])

        elif result_type == "unmatched":
            unmatched = result[result_type]
            prefix = "inter_" if region == "intergenic" else "genic_"
            for key, value in unmatched.items():
                if key.startswith(("an_" + prefix, "sn_" + prefix, "vp_" + prefix)):
                    row[key] = value

        else:
            row.update(result[result_type])

        flattened.append(row)
    return pd.DataFrame(flattened)


def process_all_files(input_path, output_path):
    """
    Process all files in the input directory and generate structured CSV outputs for all analyses.
    """
    genic_data = []
    intergenic_data = []

    loc_mapping, symbol_mapping = load_gene_mapping("C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\hgnc_complete_set.csv")

    for file_path in glob.glob(input_path):
        print(f"Processing file: {file_path}")
        results = data_process(file_path, loc_mapping, symbol_mapping)

        genic_data.append({
            "chromosome": results["chromosome"],
            "concat": results["genic_concat"],
            "concat_rate": results["genic_concat_rate"],
            "sunburst": results["genic_sunburst"],
            "unmatched": results["genic_unmatched"]
        })

        intergenic_data.append({
            "chromosome": results["chromosome"],
            "concat": results["inter_concat"],
            "concat_rate": results["inter_concat_rate"],
            "sunburst": results["inter_sunburst"],
            "unmatched": results["inter_unmatched"]
        })

    # Flatten results
    genic_concat_df = flatten_results(genic_data, "concat", "genic")
    genic_concat_rate_df = flatten_results(genic_data, "concat_rate", "genic")
    genic_sunburst_df = flatten_results(genic_data, "sunburst", "genic")

    inter_concat_df = flatten_results(intergenic_data, "concat", "intergenic")
    inter_concat_rate_df = flatten_results(intergenic_data, "concat_rate", "intergenic")
    inter_sunburst_df = flatten_results(intergenic_data, "sunburst", "intergenic")

    # Flatten unmatched results
    genic_unmatched_df = flatten_results(genic_data, "unmatched", "genic")
    inter_unmatched_df = flatten_results(intergenic_data, "unmatched", "intergenic")

    # Save to CSV
    genic_concat_df.to_csv(f"{output_path}/genic_concat_refseq.csv", index=False)
    genic_concat_rate_df.to_csv(f"{output_path}/genic_concat_rate_refseq.csv", index=False)
    genic_sunburst_df.to_csv(f"{output_path}/genic_sunburst_refseq.csv", index=False)

    inter_concat_df.to_csv(f"{output_path}/inter_concat_refseq.csv", index=False)
    inter_concat_rate_df.to_csv(f"{output_path}/inter_concat_rate_refseq.csv", index=False)
    inter_sunburst_df.to_csv(f"{output_path}/inter_sunburst_refseq.csv", index=False)

    # Save unmatched info
    genic_unmatched_df.to_csv(f"{output_path}/genic_unmatched_counts_refseq.csv", index=False)
    inter_unmatched_df.to_csv(f"{output_path}/inter_unmatched_counts_refseq.csv", index=False)


def main():
    # Define input and output paths
    input_path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\AnnoQ_data_all\\*.gz"
    #output_path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\Results_refseq_nested"
    
    #New Data
    #input_path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\NEW_ANNOQ_annotations\\chr\\*.gz"
    output_path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\Results_refseq_nested_NEW_ANNOQ"

    # Process all files and generate outputs
    process_all_files(input_path, output_path)

if __name__ == "__main__":
    main()
