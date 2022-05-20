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
        .about("Call indels and inversions from sam/bam files (sorted by chromosome/pos)")
        .arg(arg!(-i --input [FILE]              "input file name (default from 'stdin')"))
        .arg(arg!(-c --min_cov [NUMBER]              "minimum indel coverage"))
        .arg(arg!(-n --no_small_indels ... "Do not call small indels"))
        .arg(arg!(-d --debug ... "Do not call small indels"))
        .get_matches();
    // get the values
    let min_cov: isize = matches.value_of("min_cov").unwrap_or("5").parse().expect("Please give a number of minimum coverage"); // minimum coverage at a position
    eprintln!("minimum coverage is {}", min_cov);
    let in_file = matches.value_of("input").unwrap_or("stdin");
    eprintln!("Input file is {}", in_file);
    let no_small_indels = matches.is_present("no_small_indels"); // .unwrap_or(false); // default is to call
    let debug = matches.is_present("debug");

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

    let mut indel_dict = HashMap::new();
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
                if &line[0..1] != "@" {// skip headers
                    // eprintln!("{}", &line);
                    parse_line(&line, &mut indel_dict, no_small_indels, debug);
                }
                // do not accumulate data
                line.clear();
            }
            Err(err) => {
                panic!("Error - {}", err );
            }
        };
    } // end of loop
    println!("Chrom\tref_start\tref_end\talt\tmuation_size\ttype\tcoverage");
    for (key, value) in &indel_dict {
        if value >= &min_cov {
            println!("{}\t{}", &key, value);
        }
    }
}

// split cigar
fn mysplit (cigar: &str) -> (Vec<char>, Vec<isize>) {
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

fn parse_cigar (cigar: &str, ref_pos: isize, same_strand: bool, read_len: isize) -> Vec<isize> {
  let (ss1, ss2) = mysplit(cigar); //@["M", "D", "M", "S"], @["60", "5", "56", "26"]
    let mut read_pos1 = 0; // left start
    let mut read_pos2 = -1; // right end
    let ref_pos1 = ref_pos; // left start
    let mut ref_pos2 = ref_pos - 1; // right end
    let mut nmatch = 0; // number of M, if match showed up, then no more S or H
  for i in 0 .. ss1.len() {
    let num = ss2[i];
    if ss1[i] == 'M' || ss1[i] == '=' || ss1[i] == 'X'{
      read_pos2 += num;
      ref_pos2 += num;
      nmatch += 1;
    } else if ss1[i] == 'S' || ss1[i] == 'H' {
      if nmatch == 0 {
        read_pos1 += num;
        read_pos2 += num;
      }
    } else if ss1[i] == 'I'{
      read_pos2 += num;
    } else if ss1[i] == 'D' || ss1[i] == 'N' {
      ref_pos2 += num - 1;
    }
  }
  if same_strand {
    return vec![read_pos1, read_pos2, ref_pos1, ref_pos2];
  } else {
    return vec![read_len - read_pos2 - 1, read_len - read_pos1 - 1, ref_pos1, ref_pos2];
  }
}

// parse each line
fn parse_line (line: &str, map: &mut HashMap<String, isize>, no_small_indels: bool, debug: bool) {
    let ff: Vec<&str> = line.split("\t").collect();
    let flag:u32 = ff[1].parse().unwrap();
    let chrom  = ff[2];
    let pos  = ff[3].parse::<isize>().unwrap() - 1; // 0 based
    let cigar  = ff[5];
    let read_seq= ff[9];
    let strand = if (flag & 0x10) != 0 {"-"} else {"+"};
    // let read_id = ff[0];
    if cigar == "*" { // no mapping
        return;
    }
    // echo read_id
    let mut mut_size: isize;
    let mut kk: String; //String::new(); // key of indel dict
    // check small indels from mapping software
    if !no_small_indels && (cigar.contains('I') || cigar.contains('D') || cigar.contains('N')){
        if debug {
            eprintln!("Small indel Line is: {}", line);
        }
        // let all_pos1 = parse_cigar(cigar, pos, true, read_seq.len().try_into().unwrap());
      let (ss1, ss2) = mysplit(cigar); //@["M", "D", "M", "S"], @[60, 5, 56, 26]
      let mut read_pos = 0;
      let mut ref_pos = pos;
      for i in 0 .. ss1.len() {
        let num = ss2[i];
        if ss1[i] == 'M' || ss1[i] == '=' || ss1[i] == 'X' {
          read_pos += num;
          ref_pos += num;
        } else if ss1[i] == 'S' {
          read_pos += num;
        } else if ss1[i] == 'I' {
          mut_size = ss2[i];
          if read_pos == 0 || read_pos+num+1 > read_seq.len().try_into().unwrap() {
            // eprintln!("Small insertion Line is: {}", line);
            // eprintln!("read_pos, num, read_seq length is {}, {}, {}", read_pos, num, read_seq.len());
            return;
          }
          let alt_seq = &read_seq[(read_pos-1).try_into().unwrap() .. (read_pos+num+1).try_into().unwrap()];
          kk = [chrom.to_string(), ref_pos.to_string(), (ref_pos+1).to_string(), alt_seq.to_string(), mut_size.to_string(), "small_ins".to_string()].join("\t");
          let count = map.entry(kk).or_insert(0);
          *count += 1;
          read_pos += num;
        } else if ss1[i] == 'D' || ss1[i] == 'N' {
          mut_size = ss2[i];
        //   altSeq = readSeq[readPos-1 .. readPos]
          if read_pos == 0 || read_pos+1 > read_seq.len().try_into().unwrap() {
            // eprintln!("Small insertion Line is: {}", line);
            // eprintln!("read_pos, read_seq length is {}, {}", read_pos, read_seq.len());
            return;
          }
          let alt_seq = &read_seq[(read_pos-1).try_into().unwrap() .. (read_pos+1).try_into().unwrap()];
          kk = [chrom.to_string(), ref_pos.to_string(), (ref_pos+num+1).to_string(), alt_seq.to_string(), mut_size.to_string(), "small_del".to_string()].join("\t");
          let count = map.entry(kk).or_insert(0);
          *count += 1;
          ref_pos += num;
        }
      }
    }
    // check if having supplementary read
    // SA:Z:6B,185,-,54S62M26S,10,0;    
    if line.contains("SA:Z") && !cigar.contains('H'){ // only when the read has no hard clip, eg. read is raw
        let ff3: Vec<&str> = line.split("SA:Z:").collect();
        let ff2: Vec<&str> = ff3[1].split(",").collect();
        let sa_chrom  = ff2[0];
        let sa_pos  = ff2[1].parse::<isize>().unwrap() - 1; // 0-based this is close to the border on the left, may need to adjust
        let sa_strand = ff2[2];
        let sa_cigar  = ff2[3];
        let all_pos1: Vec<isize>; // = Vec::new();
        let all_pos2: Vec<isize>; // = Vec::new();
        if chrom == sa_chrom && strand == sa_strand { // potential big deletion, could be insertion too, but update later
            if sa_pos > pos { // SA is on the right
                all_pos1 = parse_cigar(cigar, pos, true, read_seq.len().try_into().unwrap());
                all_pos2 = parse_cigar(sa_cigar, sa_pos, true, read_seq.len().try_into().unwrap());
            } else {
                all_pos1 = parse_cigar(sa_cigar, sa_pos, true, read_seq.len().try_into().unwrap()); // return read_pos1, read_pos2, ref_pos1, ref_pos2
                all_pos2 = parse_cigar(cigar, pos, true, read_seq.len().try_into().unwrap());
            }
            if debug {
                eprintln!("Big deletion Line is: {}", line);
                eprintln!("Big deletions\nall_pos1 is {:?}", &all_pos1);
                eprintln!("all_pos2 is {:?}", &all_pos2);
            }
            // echo "potential big deletion"
            if all_pos1[0] > all_pos2[1] || all_pos1[3] >= all_pos2[2] || all_pos1[1] >= all_pos2[1] {return;}
            let read_pos1 = all_pos1[1];
            let ref_pos1  = all_pos1[3];
            let read_pos2 = all_pos2[0];
            let ref_pos2  = all_pos2[2];
            let shift = if read_pos1 >= read_pos2 {read_pos1 - read_pos2 + 1} else {0};
            let del_end_pos = ref_pos2 + shift;
            if ref_pos1 < del_end_pos {
                // altSeq = readSeq[readPos1 .. readPos2+shift]
                if read_pos2+shift+1 > read_seq.len().try_into().unwrap() {
                    // eprintln!("Big deletion Line is {}", line);
                    // eprintln!("all_pos1 is {:?}", &all_pos1);
                    // eprintln!("all_pos2 is {:?}", &all_pos2);
                    return;
                }
                let alt_seq = &read_seq[read_pos1.try_into().unwrap() .. (read_pos2+shift+1).try_into().unwrap()];
                mut_size = del_end_pos - ref_pos1 - 1;
                kk = [chrom.to_string(), (ref_pos1+1).to_string(), (del_end_pos+1).to_string(), alt_seq.to_string(), mut_size.to_string(), "big_del".to_string()].join("\t");
                // println!("read_id is {}", read_id);
                // println!("kk is {}", &kk);
                let count = map.entry(kk).or_insert(0);
                *count += 1;
            }
        } else if chrom == sa_chrom && strand != sa_strand {// inversions
            if sa_pos > pos { // SA is on the right
                all_pos1 = parse_cigar(cigar, pos, true, read_seq.len().try_into().unwrap()); // return read_pos1, read_pos2, ref_pos1, ref_pos2
                all_pos2 = parse_cigar(sa_cigar, sa_pos, false, read_seq.len().try_into().unwrap());
            } else {
                all_pos1 = parse_cigar(sa_cigar, sa_pos, true, read_seq.len().try_into().unwrap());
                all_pos2 = parse_cigar(cigar, pos, false, read_seq.len().try_into().unwrap());
            }
            if debug {
                eprintln!("Inversion Line is: {}", line);
                eprintln!("Inversions\nall_pos1 is {:?}", &all_pos1);
                eprintln!("Inversionss\nall_pos2 is {:?}", &all_pos2);
            }
            let mut read_pos1 = all_pos1[1];
            let mut read_pos2 = all_pos2[0];
            let mut ref_pos1  = all_pos1[3];
            let mut ref_pos2  = all_pos2[3];
            let mut shift = if read_pos1 >= read_pos2 {read_pos1 - read_pos2 + 1} else {0};
            let mut del_end_pos = ref_pos2 - shift;
            if all_pos1[0] > all_pos2[0]{ // count from right
                read_pos1 = all_pos1[0];
                read_pos2 = all_pos2[1];
                ref_pos1  = all_pos1[2];
                ref_pos2  = all_pos2[2];
                shift = if read_pos2 >= read_pos1 {read_pos2 - read_pos1 + 1} else {0};
                del_end_pos = ref_pos2 + shift;
            }
            if ref_pos1 <= del_end_pos {
                mut_size = del_end_pos - ref_pos1 - 1;
                kk = [chrom.to_string(), (ref_pos1+1).to_string(), (del_end_pos+1).to_string(), "inversion".to_string(), mut_size.to_string(), "inv".to_string()].join("\t");
                // println!("read_id is {}", read_id);
                // println!("kk is {}", &kk);
                let count = map.entry(kk).or_insert(0);
                *count += 1;
            }
        }
    }
}