import argparse
import json
from pathlib import Path
import pandas as pd

import os
from dotenv import load_dotenv
from langchain_openai.chat_models import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

def parse_fastp_report(json_file: Path) -> dict:
    """Parses a single fastp JSON report to extract a comprehensive set of QC metrics."""
    sample_name = json_file.stem
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except (IOError, json.JSONDecodeError) as e:
        print(f"Warning: Could not read or parse {json_file}. Skipping. Error: {e}")
        return {}

    # Extracting data using .get() for safety against missing keys
    summary = data.get('summary', {})
    before_filtering = summary.get('before_filtering', {})
    after_filtering = summary.get('after_filtering', {})
    adapter_cutting = data.get('adapter_cutting', {})
    duplication = data.get('duplication', {})
    insert_size = data.get('insert_size', {})

    return {
        'Sample': sample_name,
        # Before Filtering
        'Total Reads (Before)': before_filtering.get('total_reads'),
        'Total Bases (Before)': before_filtering.get('total_bases'),
        'Q20 Rate (Before)': before_filtering.get('q20_rate'),
        'Q30 Rate (Before)': before_filtering.get('q30_rate'),
        'GC Content (Before)': before_filtering.get('gc_content'),
        # After Filtering
        'Total Reads (After)': after_filtering.get('total_reads'),
        'Total Bases (After)': after_filtering.get('total_bases'),
        'Q20 Rate (After)': after_filtering.get('q20_rate'),
        'Q30 Rate (After)': after_filtering.get('q30_rate'),
        'GC Content (After)': after_filtering.get('gc_content'),
        # Other Metrics
        'Adapter Trimmed Reads': adapter_cutting.get('adapter_trimmed_reads'),
        'Duplication Rate': duplication.get('rate'),
        'Insert Size Peak': insert_size.get('peak'),
    }

def parse_featurecounts_summary(summary_file: Path) -> dict:
    """Parses the featureCounts summary file to extract alignment stats."""
    stats = {}
    try:
        with open(summary_file, 'r') as f:
            # Skip the header line
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    status, count = parts
                    stats[status] = int(count)
    except (IOError, ValueError) as e:
        print(f"Warning: Could not read or parse {summary_file}. Skipping. Error: {e}")
        return {}

    # Calculate total reads and assignment rate
    total_reads = sum(stats.values())
    assigned_reads = stats.get('Assigned', 0)
    
    if total_reads > 0:
        assignment_rate = (assigned_reads / total_reads)
    else:
        assignment_rate = 0

    stats['Total Processed Reads'] = total_reads
    stats['Assignment Rate'] = assignment_rate
    
    return stats

def generate_summary_table(fastp_data: list, featurecounts_data: dict) -> str:
    """Generates a Markdown summary from the extracted data."""
    
    # --- Fastp Summary Table ---
    fastp_summary_str = "## Quality Control Summary (fastp)\n\n"
    if fastp_data:
        # Convert list of dicts to a pandas DataFrame
        df_fastp = pd.DataFrame(fastp_data)
        # Set Sample as index
        df_fastp = df_fastp.set_index('Sample')
        # Format float values to be more readable
        for col in ['Q20 Rate (Before)', 'Q30 Rate (Before)', 'GC Content (Before)',
                    'Q20 Rate (After)', 'Q30 Rate (After)', 'GC Content (After)',
                    'Duplication Rate']:
            if col in df_fastp.columns:
                df_fastp[col] = df_fastp[col].apply(lambda x: f"{x:.2%}" if x is not None else "N/A")
        
        fastp_summary_str += df_fastp.to_markdown()
    else:
        fastp_summary_str += "No fastp data available."

    # --- FeatureCounts Summary Table ---
    fc_summary_str = "\n\n## Alignment & Counting Summary (featureCounts)\n\n"
    if featurecounts_data:
        # Create a DataFrame for featureCounts data for easy formatting
        # We can select and rename the most important metrics for the summary
        fc_metrics = {
            "Total Processed Reads": featurecounts_data.get('Total Processed Reads'),
            "Assigned Reads": featurecounts_data.get('Assigned'),
            "Assignment Rate": featurecounts_data.get('Assignment Rate'),
            "Unassigned (Multi-mapping)": featurecounts_data.get('Unassigned_MultiMapping'),
            "Unassigned (No Features)": featurecounts_data.get('Unassigned_NoFeatures'),
            "Unassigned (Ambiguity)": featurecounts_data.get('Unassigned_Ambiguity'),
        }
        df_fc = pd.DataFrame.from_dict(fc_metrics, orient='index', columns=['Value'])
        
        # Format the assignment rate as a percentage
        if 'Assignment Rate' in df_fc.index:
            rate = df_fc.loc['Assignment Rate', 'Value']
            df_fc.loc['Assignment Rate', 'Value'] = f"{rate:.2%}" if rate is not None else "N/A"

        fc_summary_str += df_fc.to_markdown()
    else:
        fc_summary_str += "No featureCounts data available."

    return fastp_summary_str + fc_summary_str

def generate_llm_report(summary_table: str, api_key: str, base_url: str, model_name: str) -> str:
    """Generates a final report using a compatible LLM, ensuring Chinese output."""
    print("--- LLM INPUT (Summary Table) ---")
    print(summary_table)
    print("---------------------------------")

    system_prompt = """你是一名生物信息学家，任务是为一份 RNA-Seq 实验撰写摘要报告。
你的任务是分析提供的数据，并生成一份全面的、人类可读的、使用【简体中文】的 Markdown 格式报告。

报告应包括以下几个部分：
1.  **总体摘要**: 对结果进行简短、高度概括的总结。
2.  **质量控制 (QC) 分析**:
    - 评论原始测序数据的质量（例如 Q20/Q30 分数、GC 含量）。
    - 通过比较“过滤前”和“过滤后”的指标，分析过滤的效果。
    - 讨论接头污染和序列重复率，并指出任何潜在的问题。
3.  **比对与计数分析**:
    - 报告整体的 reads 分配率（Assignment Rate）。
    - 讨论未分配 reads 的比例，特别关注多重比对 (multi-mapping) 和落在注释基因区域之外的 reads。高的多重比对率或“无特征”率可能预示着特定的生物学或技术现象。
4.  **潜在问题与建议**:
    - 突出显示任何质量指标较差的样本（例如，低 Q30、高重复率）。
    - 指出比对中任何潜在的问题，例如极低的分配率。
    - 如果发现任何问题，请提供简要的建议。
"""
    human_prompt = """这是提供给你的数据：
----------------
{summary_data}
----------------

请根据这些数据，用简体中文生成完整的分析报告。"""

    # Create the ChatPromptTemplate
    prompt = ChatPromptTemplate.from_messages([
        ("system", system_prompt),
        ("human", human_prompt)
    ])

    # Set up the LangChain components using the modern Chat model
    llm = ChatOpenAI(
        model_name=model_name,
        api_key=api_key,
        base_url=base_url,
        temperature=0.2
    )
    
    output_parser = StrOutputParser()

    # Build the chain using the pipe operator
    chain = prompt | llm | output_parser

    # Invoke the chain
    print("\n5. Generating final report with LLM...")
    response = chain.invoke({"summary_data": summary_table})
    
    return response

def main(fastp_dir: Path, featurecounts_dir: Path, output_file: Path):
    """
    Main function to orchestrate the parsing, summarization, and report generation.
    """
    # Explicitly load environment variables from the .env file mounted at /app/.env.
    # This is the most robust method, bypassing all environment inheritance issues.
    # The .env file is expected to be mounted as a volume during `docker run`.
    # Explicitly load environment variables from the .env file mounted at /app/.env.
    # This is the most robust method, bypassing all environment inheritance issues.
    # The .env file is expected to be mounted as a volume during `docker run`.
    load_dotenv(dotenv_path="/app/.env")
    api_key = os.getenv("OPENAI_API_KEY")
    base_url = os.getenv("OPENAI_API_BASE")
    model_name = os.getenv("OPENAI_MODEL_NAME")

    if not all([api_key, base_url, model_name]):
        raise ValueError(
            "API key, base URL, and model name are required. "
            "Please ensure OPENAI_API_KEY, OPENAI_API_BASE, and OPENAI_MODEL_NAME are set in your .env file."
        )
    print(f"1. Searching for fastp reports in: {fastp_dir}")
    fastp_files = list(fastp_dir.glob('**/*.json'))
    print(f"   Found {len(fastp_files)} reports.")

    print(f"2. Searching for featureCounts summary in: {featurecounts_dir}")
    # Assuming a single summary file, as is typical.
    featurecounts_summary_file = next(featurecounts_dir.glob('*.summary'), None)
    
    if not featurecounts_summary_file:
        raise FileNotFoundError("Could not find featureCounts summary file.")
    print(f"   Found: {featurecounts_summary_file}")

    # --- Data Parsing ---
    print("3. Parsing fastp reports...")
    all_fastp_data = [parse_fastp_report(f) for f in fastp_files if f.is_file()]
    # Filter out any empty results from failed parsing
    all_fastp_data = [d for d in all_fastp_data if d]
    
    print("--- Parsed fastp data ---")
    import pprint
    pprint.pprint(all_fastp_data)
    print("-------------------------")

    print("4. Parsing featureCounts summary...")
    featurecounts_data = parse_featurecounts_summary(featurecounts_summary_file)
    print("--- Parsed featureCounts data ---")
    pprint.pprint(featurecounts_data)
    print("-------------------------------")

    # --- Summarization ---
    summary_table = generate_summary_table(all_fastp_data, featurecounts_data)

    # --- LLM Report Generation ---
    final_report = generate_llm_report(summary_table, api_key, base_url, model_name)

    # --- Save Report ---
    output_file.write_text(final_report)
    print(f"\nSuccessfully generated report: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize RNA-Seq pipeline results using an LLM.")
    parser.add_argument("--fastp_dir", type=Path, required=True, help="Directory containing fastp JSON reports.")
    parser.add_argument("--featurecounts_dir", type=Path, required=True, help="Directory containing featureCounts output.")
    parser.add_argument("--output_file", type=Path, default="rna_seq_summary_report.md", help="Path to save the final Markdown report.")
    args = parser.parse_args()
    
    main(args.fastp_dir, args.featurecounts_dir, args.output_file)