process SEARCH_TRUNCATED_GFF {
    tag "Searching for truncated GFF in preprocessed FASTA"
    publishDir "${params.outputDir}/truncated_gff_search", mode: 'copy'
    publishDir "${projectDir}/further_analysis", pattern: "further_analysis/*", mode: 'copy'

    input:
    path truncated_gff_list
    tuple val(sampleName), path(cleaned_fasta)

    output:
    path "${sampleName}_truncated_gff_matches.txt", emit: truncated_gff_matches
    path "further_analysis/*", emit: further_analysis_files, optional: true

    script:
    """
    #!/usr/bin/env python3
    import os
    import shutil

    # Create further_analysis directory in the current work directory
    further_analysis_dir = "further_analysis"
    os.makedirs(further_analysis_dir, exist_ok=True)

    # Read the truncated GFF list
    truncated_gffs = set()
    with open("${truncated_gff_list}", 'r') as f:
        for line in f:
            gff, _ = line.strip().split('\t')
            truncated_gffs.add(gff)

    # Process the current sample
    sample_name = "${sampleName}".replace("_cleaned", "")
    matches = []

    if sample_name in truncated_gffs:
        matches.append(sample_name)
        # Copy the file to further_analysis
        new_filename = f"{sample_name}.fa"
        shutil.copy("${cleaned_fasta}", os.path.join(further_analysis_dir, new_filename))

    # Write the results
    with open("${sampleName}_truncated_gff_matches.txt", 'w') as f:
        f.write(f"Sample: {sample_name}\\n")
        f.write("Matches:\\n")
        for match in matches:
            f.write(f"{match}\\n")

    # If no matches were found, create an empty file
    if not matches:
        open(os.path.join(further_analysis_dir, "no_matches_found"), 'w').close()
    """
}