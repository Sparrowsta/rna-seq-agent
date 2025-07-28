import re

# This dictionary maps patterns in log messages to a canonical step name.
# The patterns are tried in order. The first match wins.
LOG_STEP_PATTERNS = {
    # Nextflow's own process execution messages
    re.compile(r'\[(\w+/\w+)\] process > (\w+)'): lambda m: m.group(2).upper(),
    
    # Specific tool outputs that might not be captured by Nextflow's process logger
    re.compile(r'sra-toolkit fasterq-dump'): lambda m: 'FASTERQ_DUMP',
    re.compile(r'STAR --runMode genomeGenerate'): lambda m: 'STAR_GENOME_GENERATE',
    re.compile(r'fastp'): lambda m: 'FASTP',
    re.compile(r'STAR --genomeDir'): lambda m: 'STAR_ALIGN',
    re.compile(r'featureCounts'): lambda m: 'FEATURECOUNTS',

    # Default catch-all for general Nextflow messages
    re.compile(r'\[.+] Submitted process'): lambda m: 'NEXTFLOW_SYSTEM',
}

def parse_log_line(line: str) -> tuple[str, str]:
    """
    Parses a single log line and identifies which pipeline step it belongs to.

    Args:
        line: A single line of log output.

    Returns:
        A tuple containing:
        - The identified step name (e.g., 'FASTP', 'STAR_ALIGN', 'GLOBAL').
        - The original, cleaned log line.
    """
    cleaned_line = line.strip()
    
    for pattern, step_extractor in LOG_STEP_PATTERNS.items():
        match = pattern.search(cleaned_line)
        if match:
            # If the extractor is a function, call it with the match object.
            # Otherwise, it's a simple string.
            step_name = step_extractor(match) if callable(step_extractor) else step_extractor
            return step_name, cleaned_line
            
    # If no specific pattern is matched, classify it as a general message.
    return 'GLOBAL', cleaned_line