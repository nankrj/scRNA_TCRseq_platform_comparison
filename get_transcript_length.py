import gzip
import json

def parse_gtf_for_transcripts(gtf_file, ensembl_version):
    """Parse GTF file and return transcript lengths by gene ID"""
    print(f"Parsing {ensembl_version}...")
    transcript_exons = {}
    transcript_length_dict = {}
    
    canonical_count = 0
    matched_count = 0
    
    # First pass: collect all exons
    with gzip.open(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            if fields[2] == 'exon':
                attributes = fields[8]
                transcript_id = None
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        transcript_id = attr.split('"')[1]
                        break
                
                if transcript_id:
                    exon_length = int(fields[4]) - int(fields[3]) + 1
                    if transcript_id not in transcript_exons:
                        transcript_exons[transcript_id] = 0
                    transcript_exons[transcript_id] += exon_length
    
    print(f"  Found exons for {len(transcript_exons)} transcripts")
    
    # Second pass: find canonical/basic transcripts and link to genes
    with gzip.open(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            
            if fields[2] == 'transcript':
                attributes = fields[8]
                gene_id = None
                transcript_id = None
                is_canonical = False
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                    elif attr.startswith('transcript_id'):
                        transcript_id = attr.split('"')[1]
                    elif 'tag' in attr and ('canonical' in attr or 'basic' in attr):
                        is_canonical = True
                
                if is_canonical:
                    canonical_count += 1
                    if transcript_id and transcript_id in transcript_exons:
                        transcript_length_dict[gene_id] = transcript_exons[transcript_id]
                        matched_count += 1
    
    print(f"  Canonical/basic transcripts: {canonical_count}")
    print(f"  Matched with exons: {matched_count}")
    print(f"  Final genes with lengths: {len(transcript_length_dict)}")
    
    return transcript_length_dict


# Download and parse 
print("Downloading Ensembl 109 (for Parse)...")
urllib.request.urlretrieve(
    "http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz",
    "ensembl109.gtf.gz"
)

lengths_109 = parse_gtf_for_transcripts("ensembl109.gtf.gz", "Ensembl 109")

# Save
with open('transcript_lengths_ensembl109.json', 'w') as f:
    json.dump(lengths_109, f)


print("Saved both dictionaries")