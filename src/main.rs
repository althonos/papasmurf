#![allow(unused)]

use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufRead;
use std::io::Read;
use std::io::BufReader;
use std::io::Write;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;
use std::str::FromStr;

use pssmurf::db::Database;
use pssmurf::db::DatabaseBuilder;
use pssmurf::mapper::Mapper;
use pssmurf::matrix::CsrMatrix;
use pssmurf::matrix::DokMatrix;
use pssmurf::matrix::CooMatrix;
use pssmurf::matrix::Matrix;
use pssmurf::matrix::MatrixDimensions;
use pssmurf::primer::Primer;
use pssmurf::seq::count_ambiguous;
use pssmurf::seq::dna_match;
use pssmurf::seq::mismatches;
use pssmurf::seq::reverse_complement;
use pssmurf::seq::DesambiguationIterator;
use pssmurf::utils::Interner;
use pssmurf::utils::OrderedSet;
use pssmurf::utils::Paired;
use pssmurf::io::FastqReader;
use pssmurf::io::FastaReader;
use pssmurf::utils::Rc;

use lightmotif::num::Unsigned;
use lightmotif::num::U32;
use lightmotif::pli::Maximum;
use lightmotif::pli::Encode;
use lightmotif::pli::Score;
use lightmotif::pli::Threshold;






fn main() {
    // --- FILL DATABASE

    // Create a new database builder from the given primers
    let mut builder = DatabaseBuilder::new(vec![
        Paired::new(
            Primer::new("TGGCGAACGGGTGAGTAA"),
            Primer::new(reverse_complement("CCGTGTCTCAGTCCCARTG")),
        ),
        Paired::new(
            Primer::new("ACTCCTACGGGAGGCAGC"),
            Primer::new(reverse_complement("GTATTACCGCGGCTGCTG")),
        ),
        Paired::new(
            Primer::new("GTGTAGCGGTGRAATGCG"),
            Primer::new(reverse_complement("CCCGTCAATTCMTTTGAGTT")),
        ),
        Paired::new(
            Primer::new("GGAGCATGTGGWTTAATTCGA"),
            Primer::new(reverse_complement("CGTTGCGGGACTTAACCC")),
        ),
        Paired::new(
            Primer::new("GGAGGAAGGTGGGGATGAC"),
            Primer::new(reverse_complement("AAGGCCCGGGAACGTATT")),
        ),

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
    const DB: &'static str = "gg_13_5.fasta.gz";
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
        let n_ambiguous = count_ambiguous(&read.sequence);
        if n_ambiguous == 0 {
            builder.add(&read.id, &read.sequence);
            n += 1;
        } else if n_ambiguous <= 3 {
            for dna in DesambiguationIterator::new(&read.sequence) {
                builder.add(&read.id, &dna);
            }
            n += 1;
        }
        // if i > 10000 {
        //     break
        // }
    }

    pb.finish_and_clear();
    println!("Succesfully processed {} sequences", n);

    // --- INDEX DATABASE

    println!("Building database");
    let mut db = builder.to_database();

    println!(
        "Extracted {} unique forward kmers",
        db.regions
            .iter()
            .map(|x| x.kmers.forward.columns())
            .sum::<usize>()
    );
    println!(
        "Extracted {} unique backward kmers",
        db.regions
            .iter()
            .map(|x| x.kmers.backward.columns())
            .sum::<usize>()
    );

    // --- MAP READS TO DATABASE

    // const R1: &str = "Example_L001_R1_001.fastq.gz";
    // const R2: &str = "Example_L001_R2_001.fastq.gz";
    // const R1: &str = "samples/PO49S4/PO49S4_L001_R1_001.fastq.gz";
    // const R2: &str = "samples/PO49S4/PO49S4_L001_R2_001.fastq.gz";
    // const R1: &str = "samples/MCS7/MCS7_L001_R1_001.fastq.gz";
    // const R2: &str = "samples/MCS7/MCS7_L001_R2_001.fastq.gz";
    // const R1: &str = "samples/GFS6/GFS6_L001_R1_001.fastq.gz";
    // const R2: &str = "samples/GFS6/GFS6_L001_R2_001.fastq.gz";

    const R1: &str = "raw/Q5RES023A1_20230327091114__MC_S7_R1_001.fastq";
    const R2: &str = "raw/Q5RES023A1_20230327091114__MC_S7_R2_001.fastq";

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

    let mut mapper = Mapper::new(&db);
    let mut mapped_reads = 0;

    for (i, res) in r1_reader
        .zip(r2_reader)
        .map(Paired::from)
        .enumerate()
    {
        let seq = res.map(Result::unwrap);
        if mapper.add(seq.as_ref().map(|r| r.sequence.as_str())) {
            mapped_reads += 1;
        }
        // for e in mapper.expected.iter_mut() {
        //     e.grow(1, 0);
        // }

        // let read = res.map(Result::unwrap);
        // // if read.forward.seq().len() <= builder.k || read.backward.seq().len() <= builder.k {
        // //     continue;
        // // }

        // // let mut striped = read
        // //     .as_ref()
        // //     .map(|r| pli.encode(r.sequence.as_bytes()))
        // //     .map(Result::unwrap)
        // //     .map(lightmotif::seq::EncodedSequence::<lightmotif::Dna>::from)
        // //     .map(|e| e.to_striped::<lightmotif::num::U32>());
        // // striped.as_mut().map(|s| s.configure_wrap(builder.k));

        // let seq = read.as_ref().map(|r| &r.sequence);
        // let (r, pos, primer_mismatches) = db
        // // let (r, pos) = db
        //     .regions
        //     .iter()
        //     .enumerate()
        //     // .map(|(r, region)| {
        //     //     pli.score_into(&striped.forward, &region.profile.forward, &mut scores);
        //     //     let fwd_pos = pli.argmax(&scores).unwrap();
        //     //     let fwd_score = scores[fwd_pos];

        //     //     pli.score_into(&striped.backward, &region.profile.backward, &mut scores);
        //     //     let bwd_pos = pli.argmax(&scores).unwrap();
        //     //     let bwd_score = scores[bwd_pos];

        //     //     (r, Paired::new((fwd_pos, fwd_score), (bwd_pos, bwd_score)))
        //     // })
        //     // .max_by(|(_, p1), (_, p2)| {
        //     //     (p1.forward.1 + p2.backward.1)
        //     //         .partial_cmp(&(p2.forward.1 + p2.backward.1))
        //     //         .unwrap()
        //     // })
        //     // .map(|(r, p)| (r, p.map(|x| x.0)))
        //     // .unwrap();
        //     .map(|(r, region)| {
        //         (
        //             r,
        //             (0..seq.forward.len() - region.primer.forward.len())
        //                 .map(|i| (i, region.primer.forward.mismatches(&seq.forward[i..i + region.primer.forward.len()])))
        //                 .min_by_key(|(_, s)| *s)
        //                 .unwrap(),
        //             (0..seq.backward.len() - region.primer.backward.len())
        //                 .map(|i| (i, region.primer.backward.reverse_complement().mismatches(&seq.backward[i..i + region.primer.backward.len()])))
        //                 .min_by_key(|(_, s)| *s)
        //                 .unwrap(),
        //         )
        //     })
        //     .min_by(|x, y| (x.1.1 + x.2.1).partial_cmp(&(y.1.1 + y.2.1)).unwrap())
        //     .map(|x| (x.0, Paired::new(x.1.0, x.2.0), Paired::new(x.1.1, x.2.1)))
        //     .unwrap();

        // if primer_mismatches.forward > 2 || primer_mismatches.backward > 2 {
        //     continue
        // }
        // let mut kmer = Paired::new(
        //     &seq.forward[pos.forward + db.regions[r].primer.forward.len()..],
        //     &seq.backward[pos.backward + db.regions[r].primer.backward.len()..],
        // );
        // // let mut kmer = Paired::new(
        // //     &seq.forward[pos.forward..], 
        // //     &seq.backward[pos.backward..]
        // // );
        // if kmer.forward.len() > db.k {
        //     kmer.forward = &kmer.forward[..db.k];
        // }
        // if kmer.backward.len() > db.k {
        //     kmer.backward = &kmer.backward[..db.k];
        // }
        
        // // if pos.forward + builder.k > read.forward.seq().len() {
        // //     println!("{:?}", &pos.backward);
        // //     continue;
        // // }
        // // if pos.backward + builder.k > read.backward.seq().len() {
        // //     println!("{:?}", &pos.backward);
        // //     continue;
        // // }

        // // let mut kmer = Paired::new(
        // //     if pos.forward + builder.k > read.forward.seq().len() {
        // //         &read.forward.seq()[pos.forward..]
        // //     } else {
        // //         &read.forward.seq()[pos.forward..pos.forward + builder.k]
        // //     },
        // //     if pos.backward + builder.k > read.backward.seq().len() {
        // //         &read.backward.seq()[pos.backward..]
        // //     } else {
        // //         &read.backward.seq()[pos.backward..pos.backward + builder.k]
        // //     }
        // // );

        // let mut mismatch = db.regions[r]
        //     .kmers
        //     .as_ref()
        //     .map(|block| vec![0u8; block.columns()]);
        // simd_mismatches(
        //     kmer.forward.as_bytes(),
        //     &db.regions[r].kmers.forward,
        //     &mut mismatch.forward,
        // );
        // simd_mismatches(
        //     kmer.backward.as_bytes(),
        //     &db.regions[r].kmers.backward,
        //     &mut mismatch.backward,
        // );

        // let mut mapped = false;
        // for entry in db.regions[r].entries.iter() {
        //     let mm = Paired::new(
        //         mismatch.forward[entry.kmer_index.forward],
        //         mismatch.backward[entry.kmer_index.backward],
        //     );
        //     const PE: f32 = 0.005;
        //     let ne = (mm.forward + mm.backward) as f32;
        //     let l = kmer.forward.len() + kmer.backward.len();
        //     let e = (PE / 3.0).powf(ne) * (1.0 - PE).powf(l as f32 - ne);
        //     if e > 0.0 && ne <= 2.0 {
        //         mapper.expected[r].insert(
        //             i,
        //             db.regions[r].unique_pairs[&entry.kmer_index],
        //             e,
        //         );
        //         mapped = true;
        //     }
        // }

        // if mapped {
        //     mapped_reads += 1;
        // }
        
        // if i > 1000 {
        //     break;
        // }
    }
    println!("Processed {} reads", mapper.expected[0].rows());
    pb.finish_and_clear();

    for r in 0..db.regions.len() {
        println!("[r={}] extracted: {}", r, mapper.expected[r].nnz());
    }

    // --- BUILD Q MATRIX ---
    println!("Computing Q_i,j matrix");
    let mut q_matrix = CooMatrix::<f32>::new(mapper.expected[0].rows(), db.names.len());
    for (r, region) in db.regions.iter().enumerate() {
        let e_csr = mapper.expected[r].to_csr();
        q_matrix = q_matrix + e_csr.dot(&region.matrix);
    }
    // println!("{:?}", q_matrix);


    // --- ITERATE
   
    println!("Computing π_j vector");

    let pb = indicatif::ProgressBar::new(size as u64)
        .with_style(indicatif::ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{total} ({per_sec}) {msg}")
        .unwrap());
    
    let mut pi = vec![1.0; q_matrix.columns()];
    let mut up = vec![0.0; q_matrix.columns()];
    let mut dens = vec![0.0; q_matrix.rows()];

    for it in pb.wrap_iter(0..10) {
        // println!("iteration {}", it);

        dens.fill(0.0);
        for (i, j, x) in q_matrix.iter() {
            dens[i] += x * pi[j];
        }

        up.fill(0.0);
        for (i, j, x) in q_matrix.iter() {
            if dens[i] > 0.0 {
                up[j] += *x / dens[i]
            }
        }

        for j in 0..q_matrix.columns() {
            pi[j] *= up[j] / q_matrix.rows() as f32;
        }
    }
    pb.finish_and_clear();

    println!("Computing X_j vector");
    let mut xj = vec![0.0; q_matrix.columns()];
    for j in 0..q_matrix.columns() {
        if db.amplified[j] > 0 {
            xj[j] = pi[j] / db.amplified[j] as f32;
        }
    }
    let mut tot = xj.iter().sum::<f32>();
    if tot > 0.0 {
        for j in 0..q_matrix.columns() {
            xj[j] /= tot;
        }
    }

    let reader = std::fs::File::open("gg_13_5_taxonomy.txt.gz")
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
        .collect::<HashMap<Rc<str>, Rc<str>>>();

    let mut output = std::fs::File::create("/tmp/Q5RES023A1_20230327091114__MC_S7.tsv")
        .unwrap();

    println!("Result: ({} reads)", mapped_reads);
    for j in 0..xj.len() {
        let name = &db.names[j];
        if xj[j] > 0.0 {
            writeln!(output, "{}\t{}\t{}", name, &taxonomy[&name.clone()], xj[j])
                .unwrap();
        }
        if xj[j] > 0.005 {
            println!("[{}] {}: {:?}", name, &taxonomy[&name.clone()], xj[j]);
        }
    }
}
