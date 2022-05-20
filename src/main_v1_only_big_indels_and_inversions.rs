#![allow(non_snake_case)]
use std::io::{self, BufReader, BufRead};
use std::collections::HashMap;
use clap::{arg, Command};
// use flate2::read::GzDecoder;
use std::fs::File;

fn main() {
    // parse options
    let matches = Command::new("edit_call_rust")
        .version("0.1.0")
        .author("Junli Zhang <zhjl86@gmail.com>")
        .about("Call indels and inversions from sam files")
        .arg(arg!(-i --input [FILE]              "input file name (default from 'stdin')"))
        .arg(arg!(-c --min_cov [NUMBER]              "minimum indel coverage"))
        .get_matches();
    // get the values
    let min_cov: i64 = matches.value_of("min_cov").unwrap_or("5").parse().expect("Please give a number of minimum coverage"); // minimum coverage at a position
    eprintln!("minimum coverage is {}", min_cov);
    let in_file = matches.value_of("input").unwrap_or("stdin");
    eprintln!("Input file is {}", in_file);

    let mut reader: Box<dyn BufRead> = if in_file == "stdin" {
        Box::new(BufReader::new(io::stdin()))
    // } else if in_file.ends_with(".gz") {
    //     let f = File::open(in_file).unwrap();
    //     Box::new(BufReader::new(GzDecoder::new(f)))
    } else {
        let f = File::open(in_file).unwrap();
        Box::new(BufReader::new(f))
    };

    // let mut first_line = String::new();
    // let _ = reader.read_line(&mut first_line);
    // let ncol = first_line.split("\t").count();

    let mut indelDict = HashMap::new();
    // code is from here: https://dev.to/dandyvica/different-ways-of-reading-files-in-rust-2n30
    let mut line = String::new();
    loop {
        match reader.read_line(&mut line) {
            Ok(bytes_read) => {
                // EOF: save last file address to restart from this address for next run
                if bytes_read == 0 {
                    break;
                }
                // func(&line);
                parse_line(&line, &mut indelDict);
                // do not accumulate data
                line.clear();
            }
            Err(err) => {
                panic!("Error - {}", err );
            }
        };
    } // end of loop
    println!("Chrom\tref_start\tref_end\tmuation_size\ttype\tcoverage");
    for (key, value) in &indelDict {
        if value >= &min_cov {
            println!("{}\t{}", &key, value);
        }
    }
}

// split cigar
fn mysplit (cigar: &str) -> (Vec<char>, Vec<i64>) {
    // echo cigar // 8M1D108M
    let num = "0123456789";
    let mut cc = Vec::new();
    for i in cigar.chars() {
        if !num.contains(i) {
            cc.push(i);
        }
    }
    // println!("cc is {:?}", cc);
    let mut nn = Vec::new();
    let mut tmp = String::new();
    for i in cigar.chars(){
        if num.contains(i){
            tmp.push(i);
        }
        else{
            nn.push(tmp.parse().unwrap());
            tmp = String::new();
        }
    }
    // println!("nn is {:?}", nn);
    return (cc, nn);
}

fn parseCigar (cigar: &str, refPos: i64, sameStrand: bool, readLen: i64) -> Vec<i64> {
  let (ss1, ss2) = mysplit(cigar); //@["M", "D", "M", "S"], @["60", "5", "56", "26"]
    let mut readPos1 = 0; // left start
    let mut readPos2 = -1; // right end
    let refPos1 = refPos; // left start
    let mut refPos2 = refPos - 1; // right end
    let mut nM = 0; // number of M, if match showed up, then no more S or H
  for i in 0 .. ss1.len() {
    let num = ss2[i];
    if ss1[i] == 'M' || ss1[i] == '=' || ss1[i] == 'X'{
      readPos2 += num;
      refPos2 += num;
      nM += 1;
    } else if ss1[i] == 'S' || ss1[i] == 'H' {
      if nM == 0 {
        readPos1 += num;
        readPos2 += num;
      }
    } else if ss1[i] == 'I'{
      readPos2 += num;
    } else if ss1[i] == 'D' || ss1[i] == 'N' {
      refPos2 += num - 1;
    }
  }
  if sameStrand {
    return vec![readPos1, readPos2, refPos1, refPos2];
  } else {
    return vec![readLen - readPos2 - 1, readLen - readPos1 - 1, refPos1, refPos2];
  }
}

// parse each line
fn parse_line (line: &str, map: &mut HashMap<String, i64>) {
    let ff: Vec<&str> = line.split("\t").collect();
    let flag:u32 = ff[1].parse().unwrap();
    let chrom  = ff[2];
    let pos  = ff[3].parse::<i64>().unwrap() - 1; // 0 based
    let cigar  = ff[5];
    let readSeq= ff[9];
    let strand = if (flag & 0x10) != 0 {"-"} else {"+"};
    // let readID = ff[0];
    if cigar == "*" { // no mapping
        return;
    }
    // echo readID
    let mutSize: i64;
    let kk: String; //String::new(); // key of indel dict
    // check if having supplementary read
    // SA:Z:6B,185,-,54S62M26S,10,0;    
    if line.contains("SA:Z") && !cigar.contains('H'){ // only when the read has no hard clip, eg. read is raw
        let ff3: Vec<&str> = line.split("SA:Z:").collect();
        let ff2: Vec<&str> = ff3[1].split(",").collect();
        let saChrom  = ff2[0];
        let saPos  = ff2[1].parse::<i64>().unwrap() - 1; // 0-based this is close to the border on the left, may need to adjust
        let saStrand = ff2[2];
        let saCigar  = ff2[3];
        let allPos1: Vec<i64>; // = Vec::new();
        let allPos2: Vec<i64>; // = Vec::new();
        if chrom == saChrom && strand == saStrand { // potential big deletion, could be insertion too, but update later
            if saPos > pos { // SA is on the right
                allPos1 = parseCigar(cigar, pos, true, readSeq.len().try_into().unwrap());
                allPos2 = parseCigar(saCigar, saPos, true, readSeq.len().try_into().unwrap());
            } else {
                allPos1 = parseCigar(saCigar, saPos, true, readSeq.len().try_into().unwrap()); // return readPos1, readPos2, refPos1, refPos2
                allPos2 = parseCigar(cigar, pos, true, readSeq.len().try_into().unwrap());
            }
            // echo "potential big deletion"
            if allPos1[0] > allPos2[1] || allPos1[3] >= allPos2[2] {return;} // read split fragment order reversed
            let readPos1 = allPos1[1];
            let refPos1  = allPos1[3];
            let readPos2 = allPos2[0];
            let refPos2  = allPos2[2];
            let shift = if readPos1 >= readPos2 {readPos1 - readPos2 + 1} else {0};
            let delEndPos = refPos2 + shift;
            if refPos1 <= delEndPos {
                mutSize = delEndPos - refPos1 - 1;
                kk = [chrom.to_string(), (refPos1+1).to_string(), (delEndPos+1).to_string(), mutSize.to_string(), "del".to_string()].join("\t");
                // println!("readID is {}", readID);
                // println!("kk is {}", &kk);
                let count = map.entry(kk).or_insert(0);
                *count += 1;
            }
        } else if chrom == saChrom && strand != saStrand {// inversions
            if saPos > pos { // SA is on the right
                allPos1 = parseCigar(cigar, pos, true, readSeq.len().try_into().unwrap()); // return readPos1, readPos2, refPos1, refPos2
                allPos2 = parseCigar(saCigar, saPos, false, readSeq.len().try_into().unwrap());
            } else {
                allPos1 = parseCigar(saCigar, saPos, true, readSeq.len().try_into().unwrap());
                allPos2 = parseCigar(cigar, pos, false, readSeq.len().try_into().unwrap());
            }
            let mut readPos1 = allPos1[1];
            let mut readPos2 = allPos2[0];
            let mut refPos1  = allPos1[3];
            let mut refPos2  = allPos2[3];
            let mut shift = if readPos1 >= readPos2 {readPos1 - readPos2 + 1} else {0};
            let mut delEndPos = refPos2 - shift;
            if allPos1[0] > allPos2[0]{ // count from right
                readPos1 = allPos1[0];
                readPos2 = allPos2[1];
                refPos1  = allPos1[2];
                refPos2  = allPos2[2];
                shift = if readPos2 >= readPos1 {readPos2 - readPos1 + 1} else {0};
                delEndPos = refPos2 + shift;
            }
            if refPos1 <= delEndPos {
                mutSize = delEndPos - refPos1 - 1;
                kk = [chrom.to_string(), (refPos1+1).to_string(), (delEndPos+1).to_string(), mutSize.to_string(), "inv".to_string()].join("\t");
                // println!("readID is {}", readID);
                // println!("kk is {}", &kk);
                let count = map.entry(kk).or_insert(0);
                *count += 1;
            }
        }
    }
}