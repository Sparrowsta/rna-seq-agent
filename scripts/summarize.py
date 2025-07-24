import argparse
import json
import requests
from pathlib import Path
import pandas as pd

import os
from dotenv import load_dotenv
from langchain_openai.chat_models import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

def send_webhook(url: str, report_content: str):
    """Sends the report content to the specified Feishu/Lark webhook URL as plain text."""
    if not url:
        return

    payload = {
        "msg_type": "text",
        "content": {
            "text": report_content
        }
    }

    try:
        response = requests.post(url, json=payload, headers={"Content-Type": "application/json"})
        response.raise_for_status()
        response_data = response.json()
        if response_data.get("StatusCode") == 0 or response_data.get("code") == 0:
            print(f"Successfully sent report to Feishu webhook: {url}")
        else:
            print(f"Warning: Feishu webhook returned an error. Response: {response.text}")

    except requests.exceptions.RequestException as e:
        print(f"Warning: Failed to send report to webhook {url}. Error: {e}")
    except json.JSONDecodeError:
        print(f"Warning: Could not decode JSON response from Feishu webhook. Response: {response.text}")


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

def parse_star_log(log_file: Path) -> dict:
    """Parses a STAR aligner .log file (which contains Log.final.out content)."""
    sample_name = log_file.stem
    stats = {'Sample': sample_name}
    try:
        with open(log_file, 'r') as f:
            for line in f:
                if '|' not in line:
                    continue
                # Use rsplit with maxsplit=1 to handle keys that might contain '|'
                parts = line.rsplit('|', 1)
                if len(parts) == 2:
                    key, value = [x.strip() for x in parts]
                    stats[key] = value
    except (IOError, ValueError) as e:
        print(f"Warning: Could not read or parse {log_file}. Skipping. Error: {e}")
        return {}
    return stats

def parse_featurecounts_summary(summary_file: Path) -> dict:
    """
    Parses the featureCounts summary file, which can contain one or more samples.
    It sums up the counts across all samples for a consolidated view.
    """
    try:
        # Use pandas to robustly read the tab-separated summary file
        df = pd.read_csv(summary_file, sep='\t', index_col=0)
    except (IOError, pd.errors.ParserError) as e:
        print(f"Warning: Could not read or parse {summary_file}. Skipping. Error: {e}")
        return {}

    # The actual count data is in all columns except the first (which is the index)
    # Summing across all sample columns (axis=1) for each status category
    summed_stats = df.sum(axis=1).to_dict()

    # Calculate total reads and assignment rate from the summed stats
    total_reads = sum(summed_stats.values())
    assigned_reads = summed_stats.get('Assigned', 0)
    
    if total_reads > 0:
        assignment_rate = (assigned_reads / total_reads)
    else:
        assignment_rate = 0

    # Add calculated fields to the dictionary
    summed_stats['Total Processed Reads'] = total_reads
    summed_stats['Assignment Rate'] = assignment_rate
    
    return summed_stats

def generate_summary_table(fastp_data: list, star_data: list, featurecounts_data: dict) -> str:
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

    # --- STAR Alignment Summary ---
    star_summary_str = "\n\n## Alignment Summary (STAR)\n\n"
    if star_data:
        df_star = pd.DataFrame(star_data)
        df_star = df_star.set_index('Sample')
        # Select and rename columns for clarity
        star_metrics = {
            "Number of input reads": "Input Reads",
            "Uniquely mapped reads %": "Uniquely Mapped %",
            "Number of reads mapped to multiple loci": "Multi-mapping Reads",
            "% of reads mapped to multiple loci": "Multi-mapping %",
            "Number of reads mapped to too many loci": "Too-many-loci Reads",
            "% of reads mapped to too many loci": "Too-many-loci %",
        }
        # Filter for existing columns and rename
        df_star_filtered = df_star[[col for col in star_metrics.keys() if col in df_star.columns]].rename(columns=star_metrics)
        star_summary_str += df_star_filtered.to_markdown()
    else:
        star_summary_str += "No STAR alignment data available."


    # --- FeatureCounts Summary Table ---
    fc_summary_str = "\n\n## Gene-level Counting Summary (featureCounts)\n\n"
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

    return fastp_summary_str + star_summary_str + fc_summary_str

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
    - **比对 (STAR)**: 解读 STAR 的比对统计数据。评论唯一比对率（Uniquely mapped reads %），这是成功比对的关键指标。讨论多重比对（multi-mapping）和过多位点比对（too many loci）的比例，并解释它们可能的含义（例如，重复序列、基因家族）。
    - **计数 (featureCounts)**: 报告整体的 reads 分配率（Assignment Rate）。讨论未分配 reads 的比例，并结合 STAR 的结果进行分析。例如，如果唯一比对率高但分配率低，可能意味着许多 reads 落在了基因间区。
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

def main(fastp_dir: Path, bam_dir: Path, featurecounts_dir: Path, output_file: Path, start_time: str = None, end_time: str = None, total_duration: str = None):
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

    webhook_url = os.getenv("WEBHOOK_URL")

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

    print(f"3. Searching for STAR alignment logs in: {bam_dir}")
    star_log_files = list(bam_dir.glob('**/*.log'))
    print(f"   Found {len(star_log_files)} STAR logs.")

    # --- Data Parsing ---
    print("4. Parsing fastp reports...")
    all_fastp_data = [parse_fastp_report(f) for f in fastp_files if f.is_file()]
    # Filter out any empty results from failed parsing
    all_fastp_data = [d for d in all_fastp_data if d]
    
    print("--- Parsed fastp data ---")
    import pprint
    pprint.pprint(all_fastp_data)
    print("-------------------------")

    print("5. Parsing STAR logs...")
    all_star_data = [parse_star_log(f) for f in star_log_files if f.is_file()]
    all_star_data = [d for d in all_star_data if d]
    print("--- Parsed STAR data ---")
    pprint.pprint(all_star_data)
    print("------------------------")

    print("6. Parsing featureCounts summary...")
    featurecounts_data = parse_featurecounts_summary(featurecounts_summary_file)
    print("--- Parsed featureCounts data ---")
    pprint.pprint(featurecounts_data)
    print("-------------------------------")

    # --- Summarization ---
    summary_table = generate_summary_table(all_fastp_data, all_star_data, featurecounts_data)

    # --- LLM Report Generation ---
    final_report = generate_llm_report(summary_table, api_key, base_url, model_name)

    # --- Add Time Information ---
    time_info = ""
    if start_time and end_time and total_duration:
        time_info = (
            f"# 流程运行时间\n\n"
            f"- **开始时间**: {start_time}\n"
            f"- **结束时间**: {end_time}\n"
            f"- **总耗时**: {total_duration}\n\n"
            f"---\n\n"
        )
    
    full_report = time_info + final_report

    # --- Save Report ---
    output_file.write_text(full_report)
    print(f"\nSuccessfully generated report: {output_file}")

    # --- Send Webhook ---
    if webhook_url:
        send_webhook(webhook_url, full_report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize RNA-Seq pipeline results using an LLM.")
    parser.add_argument("--fastp_dir", type=Path, required=True, help="Directory containing fastp JSON reports.")
    parser.add_argument("--bam_dir", type=Path, required=True, help="Directory containing STAR alignment logs.")
    parser.add_argument("--featurecounts_dir", type=Path, required=True, help="Directory containing featureCounts output.")
    parser.add_argument("--output_file", type=Path, default="rna_seq_summary_report.md", help="Path to save the final Markdown report.")
    parser.add_argument("--start_time", type=str, help="Pipeline start time.")
    parser.add_argument("--end_time", type=str, help="Pipeline end time.")
    parser.add_argument("--total_duration", type=str, help="Pipeline total duration.")
    args = parser.parse_args()
    
    main(args.fastp_dir, args.bam_dir, args.featurecounts_dir, args.output_file, args.start_time, args.end_time, args.total_duration)