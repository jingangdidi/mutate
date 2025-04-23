use std::path::Path;

use bio::io::fastq::{
    Reader,
    Writer,
    Record,
    //FastqRead,
};
use rand::prelude::*;

use crate::{
    error::MyError,
    utils::{my_reader, my_writer, get_output_name},
};

const ATG: &[u8;3] = &[b'A', b'T', b'G']; // ATG，加上`&`可以避免栈拷贝
const ATC: &[u8;3] = &[b'A', b'T', b'C']; // ATC，加上`&`可以避免栈拷贝
const AGC: &[u8;3] = &[b'A', b'G', b'C']; // AGC，加上`&`可以避免栈拷贝
const TGC: &[u8;3] = &[b'T', b'G', b'C']; // TGC，加上`&`可以避免栈拷贝

/// 根据指定的突变数量，对指定序列进行随机突变，将随机选取的要突变位点原始碱基，随机突变为另外3种碱基
/// 4种碱基为ATGC，例如选取的位点当前碱基是A，则从TGC中随机选取1个作为突变的新碱基
/// 没有考虑插入和缺失
fn mutate(mut sequence: Vec<u8>, ratio: f32, save_pos: bool) -> (Vec<u8>, Vec<usize>) {
    let seq_len = sequence.len();
    let num_mut = (ratio * seq_len as f32) as usize; // 计算突变位点数量
    //let mut sequence: Vec<u8> = seq.as_bytes().to_vec();
    let mut rng = rand::rng(); // 不使用随机种子
    //let mut rng = StdRng::seed_from_u64(1992); // 使用随机种子
    // 创建序列长度的连续数值索引向量，然后随机打乱
    let mut idx: Vec<usize> = (0..seq_len).collect();
    idx.shuffle(&mut rng);
    // 对打乱后的前num_mut个位置进行突变，N不算
    let mut num = 0; // 记录已突变的位点数，例如序列长100，设置突变率为100%，但是序列含有2个N，则实际突变数量为98个
    let mut mut_pos: Vec<usize> = vec![]; // 存储突变位置，1-based
    for i in 0..seq_len {
        /*
        // 获取要突变为的另外3种碱基
        let targets: Vec<u8> = NT.iter().filter(|n| n != sequence[idx[i]]).collect();
        // 随机选取一种碱基作为要突变为的碱基
        let nt = targets.choose(&mut rng).unwrap();
        // 对原始序列该位置碱基进行替换
        sequence[idx[i]] = nt;
        */

        // 随机选取一种碱基作为要突变为的碱基，并对原始序列该位置碱基进行替换
        sequence[idx[i]] = match sequence[idx[i]] {
            b'A' => {
                num += 1;
                if save_pos {
                    mut_pos.push(idx[i]+1);
                }
                *TGC.choose(&mut rng).unwrap()
            },
            b'T' => {
                num += 1;
                if save_pos {
                    mut_pos.push(idx[i]+1);
                }
                *AGC.choose(&mut rng).unwrap()
            },
            b'G' => {
                num += 1;
                if save_pos {
                    mut_pos.push(idx[i]+1);
                }
                *ATC.choose(&mut rng).unwrap()
            },
            b'C' => {
                num += 1;
                if save_pos {
                    mut_pos.push(idx[i]+1);
                }
                *ATG.choose(&mut rng).unwrap()
            },
            b'N' => b'N',
            n => panic!("sequence must be A, T, G, C, N, not {}", n as char),
        };

        // 已突变位点数达到要突变的数量则退出循环
        if num == num_mut {
            break
        }
    }
    return (sequence, mut_pos)
}

/// 对fastq文件每条read根据指定突变比例进行突变
pub fn mutate_fastq(fastq_file: &Path, outpath: &Path, ratio: f32, save_pos: bool) -> Result<(), MyError> {
    // 获取不含路径的新文件名，创建突变后的fastq文件
    let mut writer = Writer::new(my_writer(&outpath.join(get_output_name(fastq_file, ratio)))?);

    /*
    // 创建突变后的fastq文件
    // let mut writer = Writer::to_file(fastq_file.with_extension(format!("ratio_{ratio}.fastq"))).map_err(|e| MyError::WriteFastqError{file: fastq_file.to_str().unwrap().to_string(), error: e.into()})?;
    let mut writer = if fastq_file.extension().unwrap() == "gz" { // 压缩则文件名为：“原名称_ratio_xxx.gz”，例如：`test_R1.fastq.gz` -> `test_R1.fastq.ratio_0.1.gz`
        Writer::new(my_writer(&fastq_file.with_extension(format!("ratio_{ratio}.gz")))?)
    } else { // 不压缩则文件名为：“原名称_ratio_xxx.fastq”，例如：`test_R1.fastq` -> `test_R1.ratio_0.1.fastq`
        Writer::new(my_writer(&fastq_file.with_extension(format!("ratio_{ratio}.fastq")))?)
    };
    */

    /*
    // 读取原始fastq文件
    let mut reader = Reader::from_file(fastq_file).map_err(|e| MyError::ReadFastqError{file: fastq_file.to_str().unwrap().to_string(), error: e.into()})?;
    let mut record = Record::new(); // 记录每个read
    while let Ok(()) = reader.read(&mut record) {
        if record.is_empty() {
            break;
        }
        // 对该read进行突变
        let new_seq = mutate(record.seq().to_vec(), ratio);
        // 突变后的read
        let new_record = Record::with_attrs(record.id(), record.desc(), &new_seq, record.qual());
        // 保存突变后的read
        let _ = writer.write_record(&new_record);
    }
    */

    // 读取原始fastq文件
    //let mut records = Reader::from_file(fastq_file).map_err(|e| MyError::ReadFastqError{file: fastq_file.to_str().unwrap().to_string(), error: e.into()})?.records();
    let mut records = Reader::new(my_reader(fastq_file)?).records();
    while let Some(Ok(record)) = records.next() {
        // 对该read进行突变
        let (new_seq, mut mut_pos) = mutate(record.seq().to_vec(), ratio, save_pos);
        // 突变后的read
        //let new_record = Record::with_attrs(record.id(), record.desc(), &new_seq, record.qual());
        let new_record = if save_pos {
            mut_pos.sort(); // 将突变位点从小到大排序
            match record.desc() { // 该read的header的description信息
                Some(desc) => {
                    let new_desc = format!("{} {:?}", desc, mut_pos);
                    Record::with_attrs(record.id(), Some(&new_desc), &new_seq, record.qual())
                },
                None => {
                    let new_desc = format!("{:?}", mut_pos);
                    Record::with_attrs(record.id(), Some(&new_desc), &new_seq, record.qual())
                },
            }
        } else {
            Record::with_attrs(record.id(), record.desc(), &new_seq, record.qual())
        };
        // 保存突变后的read
        let _ = writer.write_record(&new_record);
    }
    Ok(())
}
