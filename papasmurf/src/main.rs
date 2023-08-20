#![allow(unused)]

use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::io::Write;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;
use std::str::FromStr;
use std::sync::Arc;

use papasmurf::Mapper;
use papasmurf::Primer;
use papasmurf::Paired;
use papasmurf::Builder;
use papasmurf::Database;
use papasmurf::io::FastaReader;
use papasmurf::io::FastqReader;
use papasmurf::seq::count_ambiguous;
use papasmurf::seq::DesambiguationIterator;

use lightmotif::num::Unsigned;
use lightmotif::num::U32;
use lightmotif::pli::Encode;
use lightmotif::pli::Maximum;
use lightmotif::pli::Score;
use lightmotif::pli::Threshold;

use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

fn main() {
    let path = std::path::PathBuf::from("/tmp/db.json");

    let db: Database = if !path.exists() {
        // --- FILL DATABASE

        // Create a new database builder from the given primers
        let mut builder = Builder::new(vec![
            Paired::new("TGGCGAACGGGTGAGTAA", "CCGTGTCTCAGTCCCARTG")
                .map(Primer::new)
                .map(Result::unwrap),
            Paired::new("ACTCCTACGGGAGGCAGC", "GTATTACCGCGGCTGCTG")
                .map(Primer::new)
                .map(Result::unwrap),
            Paired::new("GTGTAGCGGTGRAATGCG", "CCCGTCAATTCMTTTGAGTT")
                .map(Primer::new)
                .map(Result::unwrap),
            Paired::new("GGAGCATGTGGWTTAATTCGA", "CGTTGCGGGACTTAACCC")
                .map(Primer::new)
                .map(Result::unwrap),
            Paired::new("GGAGGAAGGTGGGGATGAC", "AAGGCCCGGGAACGTATT")
                .map(Primer::new)
                .map(Result::unwrap),
            // Paired::new(
            //     "TGGCGGACGGGTGAGTAA",
            //     &reverse_complement("CTGCTGCCTCCCGTAGGA"),
            // )
            // .map(Primer::new),
            // Paired::new(
            //     "TCCTACGGGAGGCAGCAG",
            //     &reverse_complement("TATTACCGCGGCTGCTGG"),
            // )
            // .map(Primer::new),
            // Paired::new(
            //     "CAGCAGCCGCGGTAATAC",
            //     &reverse_complement("CGCATTTCACCGCTACAC"),
            // )
            // .map(Primer::new),
            // Paired::new(
            //     "AGGATTAGATACCCTGGT",
            //     &reverse_complement("GAATTAAACCACATGCTC"),
            // )
            // .map(Primer::new),
            // Paired::new(
            //     "GCACAAGCGGTGGAGCAT",
            //     &reverse_complement("CGCTCGTTGCGGGACTTA"),
            // )
            // .map(Primer::new),
            // Paired::new(
            //     "AGGAAGGTGGGGATGACG",
            //     &reverse_complement("CCCGGGAACGTATTCACC"),
            // )
            // .map(Primer::new),
        ]);

        // Load reference sequences
        const DB: &'static str = "../gg_13_5.fasta.gz";
        // const DB: &'static str = "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz";
        let size = std::fs::metadata(DB).unwrap().len();
        let pb = indicatif::ProgressBar::new(size as u64)
                .with_style(indicatif::ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {bytes}/{total_bytes} ({binary_bytes_per_sec}) {msg}")
                .unwrap());
        let reader = std::fs::File::open(DB)
            .map(|r| pb.wrap_read(r))
            .map(flate2::read::GzDecoder::new)
            .map(FastaReader::from)
            .unwrap();

        // Extract reference region kmers from all sequences
        let mut n = 0;
        for (i, read) in reader.map(Result::unwrap).enumerate() {
            let seq = read.sequence.replace('U', "T");
            let n_ambiguous = count_ambiguous(&seq).unwrap();
            if n_ambiguous == 0 {
                builder.add(&read.id, &seq);
                n += 1;
            } else if n_ambiguous <= 3 {
                for dna in DesambiguationIterator::new(&seq).unwrap() {
                    builder.add(&read.id, &dna);
                }
                n += 1;
            }
        }

        pb.finish_and_clear();
        println!("Succesfully processed {} sequences", n);

        // --- INDEX DATABASE

        println!("Building database");
        let db = builder.to_database();
        let mut f = std::fs::File::create(&path).unwrap();
        serde_json::to_writer(&mut f, &db).unwrap();
        // rmp_serde::encode::write(&mut f, &db).unwrap();
        // bincode::serialize_into(f, &db).unwrap();
        db
    } else {
        println!("Loading database");
        let size = std::fs::metadata(&path).unwrap().len();
        let pb = indicatif::ProgressBar::new(size as u64)
            .with_style(indicatif::ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {bytes}/{total_bytes} ({binary_bytes_per_sec}) {msg}")
            .unwrap());
        let f = std::fs::File::open(&path)
            .map(|r| pb.wrap_read(r))
            .map(BufReader::new)
            .unwrap();
        serde_json::from_reader(f).unwrap()
        // rmp_serde::from_read(f).unwrap()
        // bincode::deserialize_from(f).unwrap()
    };

    println!(
        "Extracted {} unique forward kmers",
        db.regions
            .iter()
            .map(|x| x.unique_kmers.forward.len())
            .sum::<usize>()
    );
    println!(
        "Extracted {} unique backward kmers",
        db.regions
            .iter()
            .map(|x| x.unique_kmers.backward.len())
            .sum::<usize>()
    );

    // --- MAP READS TO DATABASE

    // const R1: &str = "Example_L001_R1_001.fastq";
    // const R2: &str = "Example_L001_R2_001.fastq";
    // const R1: &str = "samples/PO49S4/PO49S4_L001_R1_001.fastq";
    // const R2: &str = "samples/PO49S4/PO49S4_L001_R2_001.fastq";
    // const R1: &str = "samples/MCS7/MCS7_L001_R1_001.fastq";
    // const R2: &str = "samples/MCS7/MCS7_L001_R2_001.fastq";
    // const R1: &str = "samples/GFS6/GFS6_L001_R1_001.fastq";
    // const R2: &str = "samples/GFS6/GFS6_L001_R2_001.fastq";
    // const R1: &str = "raw/Q5RES023A1_20230327091114__MC_S7_R1_001.fastq";
    // const R2: &str = "raw/Q5RES023A1_20230327091114__MC_S7_R2_001.fastq";
    // const R1: &str = "samples/SPFS5/SPFS5_L001_R1_001.fastq";
    // const R2: &str = "samples/SPFS5/SPFS5_L001_R2_001.fastq";

    for (R1, R2) in [
        (
            "samples/GFS6/GFS6_L001_R1_001.fastq",
            "samples/GFS6/GFS6_L001_R2_001.fastq",
        ),
        (
            "samples/H2OS8/H2OS8_L001_R1_001.fastq",
            "samples/H2OS8/H2OS8_L001_R2_001.fastq",
        ),
        (
            "samples/MAP48S1/MAP48S1_L001_R1_001.fastq",
            "samples/MAP48S1/MAP48S1_L001_R2_001.fastq",
        ),
        (
            "samples/MAP49S3/MAP49S3_L001_R1_001.fastq",
            "samples/MAP49S3/MAP49S3_L001_R2_001.fastq",
        ),
        (
            "samples/MCS7/MCS7_L001_R1_001.fastq",
            "samples/MCS7/MCS7_L001_R2_001.fastq",
        ),
        (
            "samples/PO48S2/PO48S2_L001_R1_001.fastq",
            "samples/PO48S2/PO48S2_L001_R2_001.fastq",
        ),
        (
            "samples/PO49S4/PO49S4_L001_R1_001.fastq",
            "samples/PO49S4/PO49S4_L001_R2_001.fastq",
        ),
        (
            "samples/SPFS5/SPFS5_L001_R1_001.fastq",
            "samples/SPFS5/SPFS5_L001_R2_001.fastq",
        ),
    ]
    .iter()
    {
        println!("Creating mapper");
        let mut mapper = Mapper::new(&db)
            .with_kmer_mismatches(20)
            .with_primer_mismatches(10)
            .with_partial_hits(true);
        let mut mapped_reads = std::sync::atomic::AtomicUsize::new(0);

        let size = std::fs::metadata(R1).unwrap().len();
        let pb = indicatif::ProgressBar::new(size as u64)
            .with_style(indicatif::ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {bytes}/{total_bytes} ({binary_bytes_per_sec}) {msg}")
            .unwrap());
        let r1_reader = std::fs::File::open(R1)
            .map(|r| pb.wrap_read(r))
            // .map(flate2::read::GzDecoder::new)
            .map(FastqReader::from)
            .unwrap();
        let r2_reader = std::fs::File::open(R2)
            // .map(flate2::read::GzDecoder::new)
            .map(FastqReader::from)
            .unwrap();

        // let pli = lightmotif::Pipeline::<lightmotif::Dna, _>::avx2().unwrap();
        // let mut scores = lightmotif::pli::StripedScores::<lightmotif::num::U32>::empty();

        let reads = r1_reader
            .zip(r2_reader)
            .map(Paired::from)
            .map(|res| res.map(Result::unwrap))
            .collect::<Vec<_>>();
        pb.finish_and_clear();

        let pb = indicatif::ProgressBar::new(reads.len() as u64).with_style(
            indicatif::ProgressStyle::with_template(
                "[{elapsed_precise}] {bar:40.cyan/blue} {pos} reads/{len} reads ({per_sec}) {msg}",
            )
            .unwrap(),
        );
        reads
            .par_iter()
            .progress_with(pb)
            .enumerate()
            .for_each(|(i, read)| {
                if mapper.add(read.as_ref().map(|r| r.sequence.as_str())).unwrap() {
                    mapped_reads.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
            });

        let n_total = mapper.reads.load(std::sync::atomic::Ordering::Relaxed);
        let n_mapped = mapped_reads.load(std::sync::atomic::Ordering::Relaxed);
        println!("Processed {} reads", n_total);
        println!("Mapped {} reads", n_mapped);

        for r in 0..db.regions.len() {
            println!("[r={}] extracted: {}", r, mapper.expected[r].len());
        }

        println!("Reconstructing");
        let mut result = mapper.finish();
        println!("Refining");
        for it in 0..10 {
            result.refine();
        }

        let reader = std::fs::File::open("../gg_13_5_taxonomy.txt.gz")
            .map(flate2::read::GzDecoder::new)
            .map(std::io::BufReader::new)
            .unwrap();
        let taxonomy = reader
            .lines()
            .map(Result::unwrap)
            .map(|line| {
                let (id, lineage) = line.trim_end().split_once('\t').unwrap();
                (id.into(), lineage.into())
            })
            .collect::<HashMap<Arc<str>, Arc<str>>>();

        let p = std::path::PathBuf::from(R1);
        let p2 = p.file_name().unwrap().to_str().unwrap();
        let mut output = std::fs::File::create(&format!("/tmp/{}.tsv", p2)).unwrap();

        println!("Result:");
        writeln!(
            output,
            "#r1={} r2={} reads={} mapped={}",
            R1, R2, n_total, n_mapped
        )
        .unwrap();
        writeln!(output, "id\ttaxonomy\tselection\tproportion\tmapped").unwrap();
        for j in 0..result.x.len() {
            let name = &db.names[j];
            if result.x[j] > 0.0 {
                writeln!(
                    output,
                    "{}\t{}\t{:.9e}\t{:.9e}\t{}",
                    name,
                    &taxonomy[&name.clone()],
                    result.pi[j],
                    result.x[j],
                    result.mapped[j],
                )
                .unwrap();
            }
            if result.x[j] > 0.005 {
                println!("[{}] {}: {:?}", name, &taxonomy[&name.clone()], result.x[j]);
            }
        }
    }
}
