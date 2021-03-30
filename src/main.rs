extern crate clap;
extern crate bigtools;
extern crate bio;
extern crate rayon;
extern crate polars;


use bigtools::bigwigread::BigWigRead;
use bio::io::bed;
use clap::{App, load_yaml};
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use polars::prelude::*;
use std::fs::File;
use std::path::Path;


fn mean(numbers: &Vec<f32>) -> f32 {

    let sum: f32 = numbers.iter().sum();

    sum as f32 / numbers.len() as f32

}

fn extract_intervals(bed_fn: &str) -> Result<Vec<bed::Record>> {

    let mut bed_reader =  bed::Reader::from_file(bed_fn).unwrap();
    let mut records = Vec::new();

    for record in bed_reader.records(){

        let rec = record.expect("Error reading record");
        records.push(rec);
    }

    return Ok(records)

}

fn extract_mean_signal_for_regions(regions: &Vec<bed::Record>, bw_fn: &str) -> Result<Series> {

    let mut bw = BigWigRead::from_file_and_attach(bw_fn).expect("Can't read bigwig");
    let mut regions_means = vec![0.0; regions.len()];

    
    let mut ii = 0;
    for region in regions{

        let chrom = region.chrom();
        let start = region.start() as u32;
        let end = region.end() as u32;

        let intervals = bw.values(&chrom, start, end).unwrap();
        let interval_mean = mean(&intervals);

        regions_means[ii] = interval_mean;
        ii = ii + 1
    }

    let name = Path::new(bw_fn).file_stem().unwrap().to_str().unwrap();
    let ser = Series::new(name, regions_means);
    Ok(ser)
    
}

fn write_to_file(df: &mut DataFrame, filename: &str) -> Result<()>{

    let mut file = File::create(filename).expect("could not create file");

    CsvWriter::new(&mut file)
    .has_headers(true)
    .with_delimiter(b',')
    .finish(df)

}


fn main() {

    // Load CLI
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from(yaml).get_matches();

    // Get bed fn and bigwig fnames
    let bed_fn = matches.value_of("regions").unwrap();
    let mut bw_fnames: Vec<&str> = matches.values_of("bigwigs").unwrap().collect();

    // Series names cannot be duplicated
    bw_fnames.sort();
    bw_fnames.dedup();

    // Get output fn
    let output_fn = matches.value_of("output").unwrap();

    // Extract intervals from bed file
    println!("Extracting intervals from bed file.");
    let bed_records = extract_intervals(&bed_fn).unwrap();

    // Extract mean values from each BigWig file
    println!("Extracting signal from BigWig files.");
    let region_means_series = bw_fnames.par_iter()
                                       .progress()
                                       .map(|f| extract_mean_signal_for_regions(&bed_records, f).unwrap())
                                       .collect();
    
    // Generate a dataframe for storage
    let mut df = DataFrame::new(region_means_series).unwrap();

    // Write dataframe to csv
    println!("Writing results to {}", output_fn);
    write_to_file(&mut df, &output_fn).expect("Can't write to file");

}