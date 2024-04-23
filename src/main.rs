use std::fs::{File};
use polars::prelude::*;
use rayon::prelude::*;
use std::process::{Command, Stdio};
use std::path::{PathBuf};
use std::io::{ErrorKind,
              BufWriter,
              BufReader,
              Write,
              BufRead,
              Error};
use clap::{arg, 
           command, 
           value_parser, 
           ArgAction, 
           Command as ClapCommand };

fn cli() -> ClapCommand {
    command!()
            .arg(arg!(
                --bams <FILES> "List of bam files to convert to normalized bedgraphs."
                ).value_parser(value_parser!(String)
                ).num_args(1..)
            )
            .arg(arg!(
                --counts <FILE> "Counts of spikein genome or library to normalize with."
                ).value_parser(value_parser!(PathBuf))
            )
            .arg(arg!(
                --bigwigs <FILES> "List of paths to output bigwigs. Order MUST match that of bams."
                ).value_parser(value_parser!(String)
                ).num_args(1..)
            )
            .arg(arg!(
                --norm <FILE> "Path to write normalization files to."
                ).value_parser(value_parser!(PathBuf))
            )
}

/* Theisen method: Depreciated
fn get_scaling_factor(max_reads: f64, totalreads: f64, genome_size: f64, spikecount: f64) -> f64 {
    let max_over_total = totalreads / max_reads;
    return (max_over_total * genome_size ) / ( max_over_total * spikecount )
}
*/

fn main() -> Result<()> {
    let matches = cli().get_matches();
    let bams = matches.get_many::<String>("bams")
                      .expect("Missing bams argument")
                      .collect::<Vec<_>>();
    let bigwigs = matches.get_many::<String>("bigwigs")
                         .expect("Missing bigwigs argument")
                         .collect::<Vec<_>>();

    let df;
    let spikecount: Option<Vec<f64>>;
    {
        
        let spikeins = matches.get_one::<PathBuf>("counts").expect("Missing counts argument");
        df = CsvReader::from_path(spikeins)?
                            .infer_schema(None)
                            .has_header(false)
                            .finish().unwrap_or(DataFrame::default());
        if !df.is_empty() {                    
        spikecount = Some(df[0].f64().expect("spike count did not contain f64's")
                               .into_no_null_iter()
                               .collect());
        println!("Got spikeins, scaling bedgraphs")
        } else {
        spikecount = None;
        println!("No spikein")
        }
        
    }

    let scaling_factors: Vec<f64>;
    if let Some(spikecount_arr) = &spikecount {
    
        let min = spikecount_arr.clone()
                                 .into_iter()
                                 .reduce(f64::min)
                                 .unwrap();
        scaling_factors = spikecount_arr.into_iter()
            .map(|factor| {
                min / factor
            }).collect();

    } else {
        scaling_factors = vec![1f64; bams.len()];
    };
    println!("Got scaling factors: {:?}", scaling_factors);

    let cpus = rayon::current_num_threads() / bams.len();
 
    bams.par_iter()
            .zip(&scaling_factors)
            .zip(&bigwigs)
            .for_each(|((bam, scaling_factor), bigwig)| {


        println!("{} {} {}", bam, bigwig, scaling_factor);

        
        {
            //let bw_f = File::create(&bigwig).expect("Error parsing bedgraph file path");
            let status = Command::new("bamCoverage")
                                    .args(["-b", &bam, "--scaleFactor", &scaling_factor.to_string(), "-bs", "10", "-of", "bigwig", "-p", &format!("{cpus}"), "-o", bigwig])
                                    .status()
                                    /*.spawn().expect("Failed to spawn child")
                                    .stdout
                                    .ok_or_else(|| Error::new(ErrorKind::Other, "Could not capture stdout for bamCoverage"))*/
                                    .expect("bamCoverage init error");

            /*let mut writer = BufWriter::with_capacity(1024 * 1024 * 100, bw_f);
            let reader = BufReader::new(coverage);

            reader.lines().for_each(|l| writer.write_all((l.unwrap() + "\n").as_bytes()).expect("Error writing bigwig"));
            writer.flush().expect("Error finishing bigwig");*/
            assert!(status.success());
        }
        /*
        let status = Command::new("bedGraphToBigWig")
                                .args([&bedgraph, &chrom_info, &bigwig])
                                .status()
                                .expect("bedGraphToBigWig initialization error");
        
        */

    });
    {                            
    let norm_f = File::create(matches.get_one::<PathBuf>("norm").expect("Missing norm argument"))
                      .expect("error parsing norm_csv path");

    let mut out_buf = BufWriter::with_capacity(100 * 1024_usize.pow(2), norm_f);
    for factor in scaling_factors {
        out_buf.write_all((factor.to_string() + "\n").as_bytes()).expect("error writing norm_csv");
    }
    out_buf.flush().expect("Failed to finish writing norm_csv");
    }
    Ok(())
}