import pandas as pd
import numpy as np
import re
import glob
from collections import Counter

# === Loaders and Extractors ===

def load_data(file_path):
    columns_needed = [
        "chr", "pos", "ref", "alt", "rs_dbSNP151",
        "ANNOVAR_refseq_Gene_ID", "ANNOVAR_refseq_Closest_gene(intergenic_only)",
        "SnpEff_refseq_Gene_ID", "VEP_refseq_Gene_Name",
        "ANNOVAR_ensembl_Gene_ID", "ANNOVAR_ensembl_Closest_gene(intergenic_only)",
        "SnpEff_ensembl_Gene_ID", "VEP_ensembl_Gene_ID"
    ]
    return pd.read_csv(file_path, usecols=columns_needed, delimiter="\t", dtype=object, compression="gzip")

def extract_refseq(lists):
    processed = []
    for entry in lists:
        if not entry or entry.lower() == "nan":
            processed.append([])
            continue
        entry = re.sub(r"^CHR_START-?", "", entry)
        sub_entries = re.split(r"[|,]", entry)
        final_genes = []
        for sub in sub_entries:
            if "NONE" in sub.upper(): continue
            sub = sub.split(":")[0]
            parts = sub.split('-')
            current = parts[0]
            for i in range(1, len(parts)):
                if current[-1].isdigit() and parts[i][0].isdigit():
                    current += "-" + parts[i]
                elif re.match(r"AS\d+", parts[i]):
                    current += "-" + parts[i]
                else:
                    final_genes.append(current)
                    current = parts[i]
            final_genes.append(current)
        processed.append([g for g in final_genes if g])
    return processed

def extract_ensembl(lists):
    pattern = re.compile(r"ENSG\d+", re.IGNORECASE)
    return [pattern.findall(lst) if isinstance(lst, str) else [] for lst in lists]

# === Mapping ===

def load_gene_mapping(file_path):
    df = pd.read_csv(file_path, dtype=str).fillna("")
    return {
        "loc": {row["entrez_id"].strip(): row["uniprot_ids"].strip() for _, row in df.iterrows() if row["entrez_id"].strip()},
        "symbol": {row["symbol"].strip(): row["uniprot_ids"].strip() for _, row in df.iterrows() if row["symbol"].strip()},
        "ensembl": {row["ensembl_gene_id"].strip(): row["uniprot_ids"].strip() for _, row in df.iterrows() if row["ensembl_gene_id"].strip()}
    }

def map_refseq(genes, loc_map, sym_map):
    mapped, total_unmatched, unique_unmatched = [], 0, set()
    for gene_list in genes:
        temp, seen = [], set()
        for g in gene_list:
            mg = loc_map.get(g[3:]) if g.startswith("LOC") and g[3:].isdigit() else sym_map.get(g)
            if mg: temp.append(mg)
            else:
                if g not in seen:
                    seen.add(g)
                    total_unmatched += 1
                unique_unmatched.add(g)
        mapped.append(temp)
    return mapped, total_unmatched, len(unique_unmatched)

def map_ensembl(genes, ens_map):
    mapped, total_unmatched, unique_unmatched = [], 0, set()
    for gene_list in genes:
        temp, seen = [], set()
        for g in gene_list:
            mg = ens_map.get(g)
            if mg: temp.append(mg)
            else:
                if g not in seen:
                    seen.add(g)
                    total_unmatched += 1
                unique_unmatched.add(g)
        mapped.append(temp)
    return mapped, total_unmatched, len(unique_unmatched)

# === Sunburst Agreement ===
def sunburst_agreement(refseq, ensembl):
    result = {
        'Full Match': 0,
        'Ensembl_with_RefSeq_subset': 0,
        'RefSeq_with_Ensembl_subset': 0,
        'Partial_Agree_Overlap': 0,
        'Disjoint': 0,
        'Only RefSeq': 0,
        'Only Ensembl': 0,
        'Empty Reference': 0
    }
    for r, e in zip(refseq, ensembl):
        r, e = set(r), set(e)
        if not r and not e:
            result['Empty Reference'] += 1
        elif r == e:
            result['Full Match'] += 1
        elif r and e:
            if r < e:
                result['Ensembl_with_RefSeq_subset'] += 1
            elif e < r:
                result['RefSeq_with_Ensembl_subset'] += 1
            elif r & e:
                result['Partial_Agree_Overlap'] += 1
            else:
                result['Disjoint'] += 1
        elif r and not e:
            result['Only RefSeq'] += 1
        elif e and not r:
            result['Only Ensembl'] += 1
    return result

def agreement_rate(list1, list2):
    total_reference = 0
    agree_tool1 = 0
    agree_tool2 = 0

    for l1, l2 in zip(list1, list2):
        s1, s2 = set(l1), set(l2)
        reference = s1 | s2
        total_reference += len(reference)
        agree_tool1 += sum(1 for g in reference if g in s1)
        agree_tool2 += sum(1 for g in reference if g in s2)

    return {
        "Total Reference Size": total_reference,
        "Tool 1 Agreement Count": agree_tool1,
        "Tool 2 Agreement Count": agree_tool2
    }

# === Agreement Rate Per Tool ===
def tool_agreement_rates(ann_ref, ann_ens, snp_ref, snp_ens, vep_ref, vep_ens):
    return {
        "annovar": agreement_rate(ann_ref, ann_ens),
        "snpeff": agreement_rate(snp_ref, snp_ens),
        "vep": agreement_rate(vep_ref, vep_ens),
        "combined": agreement_rate(
            [set(a)|set(s)|set(v) for a, s, v in zip(ann_ref, snp_ref, vep_ref)],
            [set(a)|set(s)|set(v) for a, s, v in zip(ann_ens, snp_ens, vep_ens)]
        )
    }

def agreement_with_tools_sunburst(AN, SN, VP):
    """
    Tracks tool agreement for each SNP and generates labels for sunburst visualization.
    Uses the concatenated method (union of all tools) as the reference set.
    Returns a dictionary of counts.
    """
    counts = Counter()

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


# === File Processor ===
def process_file(file_path, maps):
    """
    Process one annotation file and compute agreement rates, sunburst stats, and unmatched counts
    using fallback logic across all SNPs, and stratify by Ensembl-defined region afterward.
    """
    df = load_data(file_path)
    chrom = df["chr"].iloc[0] if "chr" in df.columns else file_path.split("/")[-1]

    # Fallback: Gene ID or Closest Gene
    ann_ref = extract_refseq([
        inter if val.startswith('.') else val
        for val, inter in zip(df["ANNOVAR_refseq_Gene_ID"], df["ANNOVAR_refseq_Closest_gene(intergenic_only)"])
    ])
    ann_ens = extract_ensembl([
        inter if val.startswith('.') else val
        for val, inter in zip(df["ANNOVAR_ensembl_Gene_ID"], df["ANNOVAR_ensembl_Closest_gene(intergenic_only)"])
    ])

    snp_ref = extract_refseq(df["SnpEff_refseq_Gene_ID"])
    snp_ens = extract_ensembl(df["SnpEff_ensembl_Gene_ID"])
    vep_ref = extract_refseq(df["VEP_refseq_Gene_Name"])
    vep_ens = extract_ensembl(df["VEP_ensembl_Gene_ID"])

    # Map everything
    ann_ref_mapped, ann_ref_unmapped, ann_ref_unique = map_refseq(ann_ref, maps["loc"], maps["symbol"])
    snp_ref_mapped, snp_ref_unmapped, snp_ref_unique = map_refseq(snp_ref, maps["loc"], maps["symbol"])
    vep_ref_mapped, vep_ref_unmapped, vep_ref_unique = map_refseq(vep_ref, maps["loc"], maps["symbol"])

    ann_ens_mapped, ann_ens_unmapped, ann_ens_unique = map_ensembl(ann_ens, maps["ensembl"])
    snp_ens_mapped, snp_ens_unmapped, snp_ens_unique = map_ensembl(snp_ens, maps["ensembl"])
    vep_ens_mapped, vep_ens_unmapped, vep_ens_unique = map_ensembl(vep_ens, maps["ensembl"])

    # Classify SNPs after mapping using Ensembl
    genic_idx = ~df["ANNOVAR_ensembl_Gene_ID"].str.startswith(".", na=False)
    intergenic_idx = ~genic_idx

    # Define helper
    def split_by_mask(mask, *args):
        return tuple([ [x[i] for i in range(len(mask)) if mask[i]] for x in args ])


    # Build combined sets per tool
    ann_combined = [set(r) | set(e) for r, e in zip(ann_ref_mapped, ann_ens_mapped)]
    snp_combined = [set(r) | set(e) for r, e in zip(snp_ref_mapped, snp_ens_mapped)]
    vep_combined = [set(r) | set(e) for r, e in zip(vep_ref_mapped, vep_ens_mapped)]

    # Return results
    return {
        "chromosome": chrom,

        # Genic
        "genic_tool_rates": tool_agreement_rates(
            *split_by_mask(genic_idx, ann_ref_mapped, ann_ens_mapped, snp_ref_mapped, snp_ens_mapped, vep_ref_mapped, vep_ens_mapped)
        ),
        "genic_sunburst": sunburst_agreement(
            split_by_mask(genic_idx, [set(a) | set(s) | set(v) for a, s, v in zip(ann_ref_mapped, snp_ref_mapped, vep_ref_mapped)])[0],
            split_by_mask(genic_idx, [set(a) | set(s) | set(v) for a, s, v in zip(ann_ens_mapped, snp_ens_mapped, vep_ens_mapped)])[0]
        ),
        "genic_tool_sunburst": agreement_with_tools_sunburst(
            split_by_mask(genic_idx, ann_combined)[0],
            split_by_mask(genic_idx, snp_combined)[0],
            split_by_mask(genic_idx, vep_combined)[0]
        ),
        "genic_unmatched": {
            "an_ref_unmapped": ann_ref_unmapped, "an_ref_unique": ann_ref_unique,
            "sn_ref_unmapped": snp_ref_unmapped, "sn_ref_unique": snp_ref_unique,
            "vp_ref_unmapped": vep_ref_unmapped, "vp_ref_unique": vep_ref_unique,
            "an_ens_unmapped": ann_ens_unmapped, "an_ens_unique": ann_ens_unique,
            "sn_ens_unmapped": snp_ens_unmapped, "sn_ens_unique": snp_ens_unique,
            "vp_ens_unmapped": vep_ens_unmapped, "vp_ens_unique": vep_ens_unique
        },

        # Intergenic
        "inter_tool_rates": tool_agreement_rates(
            *split_by_mask(intergenic_idx, ann_ref_mapped, ann_ens_mapped, snp_ref_mapped, snp_ens_mapped, vep_ref_mapped, vep_ens_mapped)
        ),
        "inter_sunburst": sunburst_agreement(
            split_by_mask(intergenic_idx, [set(a) | set(s) | set(v) for a, s, v in zip(ann_ref_mapped, snp_ref_mapped, vep_ref_mapped)])[0],
            split_by_mask(intergenic_idx, [set(a) | set(s) | set(v) for a, s, v in zip(ann_ens_mapped, snp_ens_mapped, vep_ens_mapped)])[0]
        ),
        "inter_tool_sunburst": agreement_with_tools_sunburst(
            split_by_mask(intergenic_idx, ann_combined)[0],
            split_by_mask(intergenic_idx, snp_combined)[0],
            split_by_mask(intergenic_idx, vep_combined)[0]
        ),
        "inter_unmatched": {
            "an_ref_unmapped": ann_ref_unmapped, "an_ref_unique": ann_ref_unique,
            "sn_ref_unmapped": snp_ref_unmapped, "sn_ref_unique": snp_ref_unique,
            "vp_ref_unmapped": vep_ref_unmapped, "vp_ref_unique": vep_ref_unique,
            "an_ens_unmapped": ann_ens_unmapped, "an_ens_unique": ann_ens_unique,
            "sn_ens_unmapped": snp_ens_unmapped, "sn_ens_unique": snp_ens_unique,
            "vp_ens_unmapped": vep_ens_unmapped, "vp_ens_unique": vep_ens_unique
        }
    }



def flatten_results(data, result_type, region):
    flattened = []
    for result in data:
        row = {"chromosome": result["chromosome"], "region": region}

        if result_type in ["sunburst", "tool_sunburst"]:
            row.update(result[result_type])

        elif result_type == "unmatched":
            unmatched = result[result_type]
            for key, value in unmatched.items():
                row[key] = value


        elif result_type == "tool_rates":
            for tool_name, metrics in result[result_type].items():
                row[f"{tool_name}_total_reference"] = metrics["Total Reference Size"]
                row[f"{tool_name}_refseq_agree"] = metrics["Tool 1 Agreement Count"]
                row[f"{tool_name}_ensembl_agree"] = metrics["Tool 2 Agreement Count"]

        else:
            row.update(result[result_type])

        flattened.append(row)
    return pd.DataFrame(flattened)


# === Flatten Results ===
def process_all_files(input_path, output_path):
    """
    Process all files in the input directory and generate structured CSV outputs for all analyses.
    """
    genic_data = []
    intergenic_data = []

    gene_maps = load_gene_mapping("C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\hgnc_complete_set.csv")

    for file_path in glob.glob(input_path):
        print(f"Processing file: {file_path}")
        results = process_file(file_path, gene_maps)

        genic_data.append({
            "chromosome": results["chromosome"],
            "tool_rates": results["genic_tool_rates"],
            "sunburst": results["genic_sunburst"],
            "tool_sunburst": results["genic_tool_sunburst"],
            "unmatched": results["genic_unmatched"]
        })

        intergenic_data.append({
            "chromosome": results["chromosome"],
            "tool_rates": results["inter_tool_rates"],
            "sunburst": results["inter_sunburst"],
            "tool_sunburst": results["inter_tool_sunburst"],
            "unmatched": results["inter_unmatched"]
        })

    # Flatten results
    genic_sunburst_df = flatten_results(genic_data, "sunburst", "genic")
    inter_sunburst_df = flatten_results(intergenic_data, "sunburst", "intergenic")

    genic_tool_df = flatten_results(genic_data, "tool_rates", "genic")
    inter_tool_df = flatten_results(intergenic_data, "tool_rates", "intergenic")

    genic_unmatched_df = flatten_results(genic_data, "unmatched", "genic")
    inter_unmatched_df = flatten_results(intergenic_data, "unmatched", "intergenic")

    # Save to CSV
    genic_sunburst_df.to_csv(f"{output_path}/genic_sunburst_r_e.csv", index=False)
    inter_sunburst_df.to_csv(f"{output_path}/inter_sunburst_r_e.csv", index=False)

    genic_tool_df.to_csv(f"{output_path}/genic_tool_rates_r_e.csv", index=False)
    inter_tool_df.to_csv(f"{output_path}/inter_tool_rates_r_e.csv", index=False)

    genic_unmatched_df.to_csv(f"{output_path}/genic_unmatched_counts_r_e.csv", index=False)
    inter_unmatched_df.to_csv(f"{output_path}/inter_unmatched_counts_r_e.csv", index=False)
    
    genic_tool_sunburst_df = flatten_results(genic_data, "tool_sunburst", "genic")
    inter_tool_sunburst_df = flatten_results(intergenic_data, "tool_sunburst", "intergenic")
    
    genic_tool_sunburst_df.to_csv(f"{output_path}/genic_tool_sunburst_r_e.csv", index=False)
    inter_tool_sunburst_df.to_csv(f"{output_path}/inter_tool_sunburst_r_e.csv", index=False)


# === Main Entrypoint ===
def main():
    #Old Data
    input_path = "C:/Users/bryan/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/AnnoQ_data_all/*.gz"
    output_path = "C:/Users/bryan/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/UniProt-Results_Ensembl_RefSeq_with_Ensembl_regions"

    
    process_all_files(input_path, output_path)

if __name__ == "__main__":
    main()