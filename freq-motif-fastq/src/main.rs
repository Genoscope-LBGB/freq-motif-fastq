use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::Command;
use uuid::Uuid;

/// Analyze FASTQ files (including gzip) and generate statistics on motifs and low-complexity bases
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTQ file (supports gzip)
    #[arg(short, long)]
    input: String,

    /// Output directory to save results (default: unique directory in current directory)
    #[arg(short, long)]
    output_dir: Option<String>,

    /// Maximum number of reads to analyze
    #[arg(short, long, default_value_t = 100_000)]
    max_reads: usize,

    /// Minimum proportion to consider (in percentage)
    #[arg(short, long, default_value_t = 50.0)]
    ratio: f64,

    /// Number of initial reads to skip
    #[arg(short = 'S', long, default_value_t = 10_000)]
    skip: usize,
}

fn main() {
    let args = Args::parse();

    let input_file = args.input;
    let max_reads = args.max_reads;
    let min_proportion = args.ratio / 100.0;
    let skip_reads = args.skip;

    // Determine output directory
    let output_dir = match args.output_dir {
        Some(dir) => PathBuf::from(dir),
        None => {
            let current_dir = env::current_dir().expect("Failed to get current directory");
            current_dir.join(format!("freq_motif_{}", Uuid::new_v4()))
        }
    };

    // Create the output directory if it doesn't exist
    fs::create_dir_all(&output_dir).expect("Failed to create output directory");

    let output_csv = output_dir.join("freq-motif.csv");
    let output_graph = output_dir.join("barplot_freq-motif.png");
    let fasta_file = output_dir.join("temp.fasta");

    eprintln!("Opening the input file: {}", input_file);
    let reader = open_fastq(&input_file).expect("Error opening the FASTQ file");

    eprintln!("Skipping the first {} reads...", skip_reads);
    let (motif_counts, total_reads, read_lengths) =
        process_reads_and_write_fasta(reader, &fasta_file, max_reads, min_proportion, skip_reads);

    eprintln!("Running SDUST...");
    let status = Command::new("sdust")
        .arg(&fasta_file)
        .output()
        .expect("Failed to execute SDUST");

    if !status.status.success() {
        eprintln!("SDUST execution failed.");
        std::process::exit(1);
    }

    let dust_output_data = String::from_utf8(status.stdout).expect("Failed to read SDUST output");
    let low_complexity_reads = parse_dust_output(&dust_output_data, &read_lengths, min_proportion);

    eprintln!("Sorting results...");

    let mut proportions: HashMap<String, f64> = initialize_all_motifs();

    for (motif, count) in motif_counts {
        proportions.insert(motif.clone(), (count as f64 / total_reads as f64) * 100.0);
    }

    // Add low-complexity result
    let low_complexity_proportion = (low_complexity_reads as f64 / total_reads as f64) * 100.0;
    proportions.insert("LowComplexity".to_string(), low_complexity_proportion);

    // Sort by proportion descending
    let mut sorted_proportions: Vec<(String, f64)> = proportions.into_iter().collect();
    sorted_proportions.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    eprintln!("Saving results to CSV: {}", output_csv.display());
    save_to_csv(&output_csv, &sorted_proportions).expect("Error saving the CSV file");

    eprintln!("Running R script for barplot...");
    let status = Command::new("generate_barplot.R")
        .arg(output_csv.as_os_str())
        .arg(output_graph.as_os_str())
        .arg(args.ratio.to_string()) // Pass the -r parameter
        .status()
        .expect("Failed to run Rscript");

    if !status.success() {
        eprintln!("Rscript execution failed");
        std::process::exit(1);
    }

    eprintln!("Cleaning up temporary files...");
    fs::remove_file(&fasta_file).expect("Failed to remove temporary FASTA file");

    eprintln!(
        "Analysis completed successfully. Results saved to '{}' and '{}'.",
        output_csv.display(),
        output_graph.display()
    );
}

/// Initializes all possible motifs with 0 proportions
fn initialize_all_motifs() -> HashMap<String, f64> {
    let mut all_motifs = HashMap::new();
    let bases = ["A", "T", "G", "C"];

    for &b1 in &bases {
        for &b2 in &bases {
            all_motifs.insert(format!("{}{}", b1, b2), 0.0);
        }
    }

    for &b1 in &bases {
        for &b2 in &bases {
            for &b3 in &bases {
                all_motifs.insert(format!("{}{}{}", b1, b2, b3), 0.0);
            }
        }
    }

    all_motifs
}

/// Opens a FASTQ file, handling both plain and gzipped formats
fn open_fastq(filename: &str) -> io::Result<Box<dyn BufRead>> {
    if filename.ends_with(".gz") {
        let file = File::open(filename)?;
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        let file = File::open(filename)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Processes reads, writes to a temporary FASTA file, and counts motifs
fn process_reads_and_write_fasta<R: BufRead>(
    reader: R,
    fasta_file: &PathBuf,
    max_reads: usize,
    min_proportion: f64,
    skip_reads: usize,
) -> (HashMap<String, u64>, usize, HashMap<String, usize>) {
    let mut motif_counts: HashMap<String, u64> = HashMap::new();
    let mut total_reads = 0;
    let mut read_lengths: HashMap<String, usize> = HashMap::new();
    let mut fasta_writer = File::create(fasta_file).expect("Failed to create temporary FASTA file");

    let mut lines = reader.lines();

    // Skip the first `skip_reads` reads
    for _ in 0..skip_reads {
        lines.next(); // Header
        lines.next(); // Sequence
        lines.next(); // '+'
        lines.next(); // Quality
    }

    while let Some(Ok(header)) = lines.next() {
        if let Some(Ok(sequence)) = lines.next() {
            if total_reads >= max_reads {
                break;
            }

            let length = sequence.len();
            if length < 2 {
                continue;
            }

            // Increment total_reads only for valid sequences
            total_reads += 1;

            if total_reads % 10_000 == 0 {
                eprintln!("Processed {} reads...", total_reads);
            }

            let read_id = &header[1..];

            // Write to FASTA for SDUST
            writeln!(fasta_writer, ">{}", read_id).expect("Failed to write FASTA header");
            writeln!(fasta_writer, "{}", sequence).expect("Failed to write FASTA sequence");
            read_lengths.insert(read_id.to_string(), length);

            let mut motif_frequencies: HashMap<String, u64> = HashMap::new();

            for i in 0..length {
                if i + 2 <= length {
                    let dinucleotide = &sequence[i..i + 2];
                    *motif_frequencies
                        .entry(dinucleotide.to_string())
                        .or_insert(0) += 1;
                }

                if i + 3 <= length {
                    let trinucleotide = &sequence[i..i + 3];
                    *motif_frequencies
                        .entry(trinucleotide.to_string())
                        .or_insert(0) += 1;
                }
            }

            let total_dinucleotides = (length - 1) as f64;
            let total_trinucleotides = (length - 2) as f64;

            for (motif, count) in motif_frequencies {
                let proportion = if motif.len() == 2 {
                    (count as f64) / total_dinucleotides
                } else {
                    (count as f64) / total_trinucleotides
                };

                if proportion > min_proportion {
                    *motif_counts.entry(motif).or_insert(0) += 1;
                }
            }
        }

        lines.next(); // Skip '+'
        lines.next(); // Skip quality
    }

    eprintln!("Total reads processed: {}", total_reads);
    (motif_counts, total_reads, read_lengths)
}

/// Parses the SDUST output to calculate low-complexity reads
fn parse_dust_output(
    dust_data: &str,
    read_lengths: &HashMap<String, usize>,
    min_proportion: f64,
) -> usize {
    let mut low_complexity_reads = 0;
    let mut masked_bases: HashMap<String, usize> = HashMap::new();

    for line in dust_data.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 3 {
            let read_name = parts[0].to_string();
            let start: usize = parts[1].parse().expect("Invalid start position");
            let end: usize = parts[2].parse().expect("Invalid end position");

            let length = end - start + 1;
            *masked_bases.entry(read_name).or_insert(0) += length;
        }
    }

    for (read_name, masked) in masked_bases {
        if let Some(&total_length) = read_lengths.get(&read_name) {
            let proportion = (masked as f64) / (total_length as f64);
            if proportion > min_proportion {
                low_complexity_reads += 1;
            }
        }
    }

    low_complexity_reads
}

/// Saves the results to a CSV file
fn save_to_csv(
    output_file: &PathBuf,
    data: &[(String, f64)],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(output_file)?;
    writeln!(file, "Motif,Proportion")?;
    for (motif, proportion) in data {
        writeln!(file, "{},{:.4}", motif, proportion)?;
    }
    Ok(())
}
