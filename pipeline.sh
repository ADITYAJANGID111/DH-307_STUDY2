#!/bin/bash

SECONDS=0

# Directory where your data is located
data_dir="/home/user/Desktop/Aditya_DH307/study2/data"

# Directory for FastQC output
fastqc_output_dir="/home/user/Desktop/Aditya_DH307/study2/fastqc_output"

# Directory for trimmed data
trimmed_output_dir="/home/user/Desktop/Aditya_DH307/study2/trimmed_data"

# Directory for output HISAT2
hisat2_output_dir="/home/user/Desktop/Aditya_DH307/study2/HISAT2_output"

# Directory for StringTie output
stringtie_output_dir="/home/user/Desktop/Aditya_DH307/study2/StringTie_output"

# Reference annotation file
reference_annotation="/home/user/Desktop/Aditya_DH307/study2/grch38/genome.gtf"
# Step 1: Run FastQC on data in each subfolder
for sample_dir in "$data_dir"/*; do
    sample_name=$(basename "$sample_dir")
    
    # Create the output directory if it doesn't exist
    mkdir -p "$fastqc_output_dir/$sample_name"
    
    fastqc --threads 16 "$sample_dir"/* -o "$fastqc_output_dir/$sample_name"
done

# Step 2: Trim data using Trimmomatic
for sample_dir in "$data_dir"/*; do
    sample_name=$(basename "$sample_dir")
    
    # Create a subdirectory for each sample in the trimmed output directory
    mkdir -p "$trimmed_output_dir/$sample_name"

    for file1 in "$sample_dir"/*_1.fastq.gz; do
        file2="${file1%_1.fastq.gz}_2.fastq.gz"
        base_name=$(basename "$file1" _1.fastq.gz)

        # Run Trimmomatic for each pair of input files
        java -jar /home/user/Desktop/jagmeet_dh307/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 "$file1" "$file2" \
            "$trimmed_output_dir/$sample_name/${base_name}_output_1_paired.fastq.gz" "$trimmed_output_dir/$sample_name/${base_name}_output_1_unpaired.fastq.gz" \
            "$trimmed_output_dir/$sample_name/${base_name}_output_2_paired.fastq.gz" "$trimmed_output_dir/$sample_name/${base_name}_output_2_unpaired.fastq.gz" \
            TRAILING:10 -phred33
        
        rm "$trimmed_output_dir/$sample_name/${base_name}_output_1_unpaired.fastq.gz"
        rm "$trimmed_output_dir/$sample_name/${base_name}_output_2_unpaired.fastq.gz"
    done
done

# Step 3: Run HISAT2 and create .bam files sample-wise
for sample_dir in "$trimmed_output_dir"/*; do
    sample_name=$(basename "$sample_dir")

    # Create a subdirectory for each sample in the HISAT2 output directory
    mkdir -p "$hisat2_output_dir/$sample_name"

    for file1 in "$sample_dir"/*_1_paired.fastq.gz; do
        file2="${file1%_1_paired.fastq.gz}_2_paired.fastq.gz"
        base_name=$(basename "$file1" _1_paired.fastq.gz)

        # Run HISAT2 for each pair of trimmed FASTQ files and create a .bam file
        hisat2 -q --threads 16 -x grch38/genome -1 "$file1" -2 "$file2" -S "$hisat2_output_dir/$sample_name/${base_name}.sam"
        samtools view -bS "$hisat2_output_dir/$sample_name/${base_name}.sam" -o "$hisat2_output_dir/$sample_name/${base_name}.bam"
        
        rm "$trimmed_output_dir/$sample_name/${base_name}_1_paired.fastq.gz"
        rm "$trimmed_output_dir/$sample_name/${base_name}_2_paired.fastq.gz"
        rm "$hisat2_output_dir/$sample_name/${base_name}.sam"
        
        # Sort the .bam file by position
        samtools sort --threads 16 "$hisat2_output_dir/$sample_name/${base_name}.bam" -o "$hisat2_output_dir/$sample_name/${base_name}.sorted.bam"
        mv "$hisat2_output_dir/$sample_name/${base_name}.sorted.bam" "$hisat2_output_dir/$sample_name/final_merged.bam"
        samtools index "$hisat2_output_dir/$sample_name/final_merged.bam"
    done
done

echo "HISAT2 finished processing all input files."

# Step 4: Run StringTie for transcript assembly
for sample_dir in "$hisat2_output_dir"/*; do
    # Check if the directory contains the final merged BAM file
    if [ -f "$sample_dir/final_merged.bam" ]; then
        sample_name=$(basename "$sample_dir")

        # Create a subdirectory for each sample in the StringTie output directory
        mkdir -p "$stringtie_output_dir/$sample_name"

        # Debugging: Print information about the files being processed
        echo "Processing sample: $sample_name"

        # Run StringTie for the merged BAM file to generate GTF file
        stringtie "$sample_dir/final_merged.bam" -o "$stringtie_output_dir/$sample_name/${sample_name}_transcripts.gtf" -G "$reference_annotation" -e

        # Debugging: Print the contents of the output directory
        echo "Contents of $stringtie_output_dir/$sample_name:"
        ls -l "$stringtie_output_dir/$sample_name"
    else
        echo "Skipping directory $sample_dir as it does not contain the final merged BAM file."
    fi
done
