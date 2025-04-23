use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::{
    Compression, // 压缩级别
    read::GzDecoder,
    write::GzEncoder,
};

use crate::error::MyError;

/// 读取fastq或fastq.gz，这样可以返回fastq或fastq.gz的reader
/// 参考：https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561/2
/// 参考：https://github.com/rust-lang/flate2-rs/issues/393
pub fn my_reader(file: &Path) -> Result<Box<dyn BufRead>, MyError> {
    let opened_file = File::open(&file).map_err(|e| MyError::ReadFileError{file: file.to_str().unwrap().to_string(), error: e})?;
    if file.extension().unwrap() == "gz" {
        //Box::new(BufReader::with_capacity(8 * 1024, GzDecoder::new(opened_file)))
        Ok(Box::new(BufReader::new(GzDecoder::new(opened_file))))
    } else {
        //Box::new(BufReader::with_capacity(8 * 1024, opened_file))
        Ok(Box::new(BufReader::new(opened_file)))
    }
}

/// 保存fastq或fastq.gz文件
/// 参考：https://users.rust-lang.org/t/write-to-normal-or-gzip-file-transparently/35561/2
pub fn my_writer(file: &Path) -> Result<Box<dyn Write>, MyError> {
    let created_file = File::create(file).map_err(|e| MyError::CreateFileError{file: file.to_str().unwrap().to_string(), error: e})?;
    if file.extension().unwrap() == "gz" {
        Ok(Box::new(BufWriter::new(GzEncoder::new(created_file, Compression::default()))))
    } else {
        //Ok(Box::new(BufWriter::with_capacity(128 * 1024, created_file)))
        Ok(Box::new(BufWriter::new(created_file)))
    }
}

/// 计算耗时，并返回字符串
pub fn my_time(info: &str, t: u128) -> String {
    if t / 600000000000 > 0 { // 大于10分钟则不显示秒
        format!("{}: {} (min)", info, t / 60000000000)
    } else if t / 1000000000 > 0 {
        format!("{}: {} (s)", info, t / 1000000000)
    } else if t / 1000000 > 0 {
        format!("{}: {} (ms)", info, t / 1000000)
    } else if t / 1000 > 0 {
        format!("{}: {} (us)", info, t / 1000)
    } else {
        format!("{}: {} (ns)", info, t)
    }
}

/// 根据输入的thread，将待处理数量n分为指定thread个子部分，返回二维数组，存储每组的起始和终止index
/// 若指定n小于thread，则每个一组
/// 若总数n/thread有余数，则将余数分给每个组1个
pub fn split_2_group(thread: usize, n: usize) -> Vec<Vec<usize>> {
    let mut group: Vec<Vec<usize>> = Vec::new();
    if n < thread {
        for i in 0..n {
            group.push(vec![i, i+1]);
        }
    } else {
        let tmp_num1 = n/thread; // 每组内数量
        let mut tmp_num2 = n%thread; // 若余数>0，则将余数平摊给每个组1个，此时每组内数量分为2种，数量差1个
        let mut m = 0;
        for i in 0..thread {
            if tmp_num2 > 0 { // 这几个组加1，平摊余数
                group.push(vec![i*tmp_num1+m, (i+1)*tmp_num1+m+1]);
                tmp_num2 -= 1;
                m = m+1;
            } else {
                group.push(vec![i*tmp_num1+m, (i+1)*tmp_num1+m]);
            }
        }
    }
    //println!("\tsplit matrix to {} part(s): {:?}", thread, group);
    group
}

/// 可能的格式后缀，前12种判断是否是双端fastq，后6种判断是否是单端fastq
const POSSIBLE_SUFFIX: &[&str;18] = &[
    "_R1.fastq", "_R1.fastq.gz",
    "_R1.fq",    "_R1.fq.gz",
    "_R1.txt",   "_R1.txt.gz",
    "_R2.fastq", "_R2.fastq.gz",
    "_R2.fq",    "_R2.fq.gz",
    "_R2.txt",   "_R2.txt.gz",
    ".fastq",    ".fastq.gz",
    ".fq",       ".fq.gz",
    ".txt",      ".txt.gz",
];

/// 从fastq文件名中提取样本名，插入“_ratio_xxx”，并保持格式后缀不变，返回不含路径的新文件名
pub fn get_output_name(fastq_file: &Path, ratio: f32) -> String {
    let tmp_name_with_siffux = fastq_file.file_name().unwrap().to_str().unwrap();
    for n in POSSIBLE_SUFFIX {
        if tmp_name_with_siffux.ends_with(n) {
            // test_R1.fastq.gz --> test_ratio_0.1_R1.fastq.gz
            // test.fastq.gz    --> test_ratio_0.1.fastq.gz
            return format!("{}_ratio_{}{}", tmp_name_with_siffux.strip_suffix(n).unwrap(), ratio, n)
        }
    }
    tmp_name_with_siffux.to_string()
}
