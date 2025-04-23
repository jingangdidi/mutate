use std::fs::create_dir_all;
use std::path::PathBuf;

use argh::FromArgs;

/// error: 定义的错误类型，用于错误传递
use crate::{
    error::MyError,
};

//----------------------------------------------------------------------------------------------------------------
#[derive(FromArgs)]
/// mutate fastq
struct Paras {
    /// fastq file (support gz compression), separate with commas
    #[argh(option, short = 'f')]
    fastq: String,

    /// mutation ratio, 0~1
    #[argh(option, short = 'r')]
    ratio: f32,

    /// max thread number, default: 4
    #[argh(option, short = 't')]
    thread: Option<usize>,

    /// add mutation position in header description
    #[argh(switch, short = 'p')]
    pos: bool,

    /// output path, default: ./
    #[argh(option, short = 'o')]
    outpath: Option<String>,
}

//----------------------------------------------------------------------------------------------------------------
/// 存储解析后的命令行参数
///#[derive(Debug, Default)]
pub struct ParsedParas {
    pub fastq:   Vec<PathBuf>, // fastq文件，多个之间逗号间隔，支持gz压缩
    pub ratio:   f32,          // 突变比例，0~1
    pub pos:     bool,         // 指定是否将突变位点的位置保存在read的header中，例如：`@LH00128:306:22MTVJLT4:2:1101:25772:1056 1:N:0:CAAGGTGA+CGGCTATG` -> `@LH00128:306:22MTVJLT4:2:1101:25772:1056 1:N:0:CAAGGTGA+CGGCTATG [1, 3, 6, 8, 9]`
    pub thread:  usize,        // 线程数，默认4
    pub outpath: PathBuf,      // 输出结果路径，不存在则自动创建，默认./
}

//----------------------------------------------------------------------------------------------------------------
/// 解析参数
pub fn parse_para() -> Result<ParsedParas, MyError> {
    let para: Paras = argh::from_env();
    let out: ParsedParas = ParsedParas{
        fastq: { // fastq文件，多个之间逗号间隔，支持gz压缩
            let mut fastq_files: Vec<PathBuf> = vec![];
            for f in para.fastq.split(",") {
                let tmp_file = PathBuf::from(f);
                if !(tmp_file.exists() && tmp_file.is_file()) {
                    return Err(MyError::FileNotExistError{file: f.to_string()})
                }
                fastq_files.push(tmp_file);
            }
            fastq_files
        },
        ratio: { // 突变比例，0~1
            if para.ratio <= 0.0 {
                return Err(MyError::ParaError{para: format!("-r mutation ratio must > 0, not {}", para.ratio)})
            } else if para.ratio > 1.0 {
                return Err(MyError::ParaError{para: format!("-r mutation ratio must <= 1, not {}", para.ratio)})
            }
            para.ratio
        },
        pos: para.pos, // 指定是否将突变位点的位置（逗号间隔）保存在read的header中，例如：`@LH00128:306:22MTVJLT4:2:1101:25772:1056 1:N:0:CAAGGTGA+CGGCTATG` -> `@LH00128:306:22MTVJLT4:2:1101:25772:1056 1:N:0:CAAGGTGA+CGGCTATG 1,3,6,8,9`
        thread: match para.thread { // 线程数，默认4
            Some(thread) => {
                if thread == 0 {
                    return Err(MyError::ParaError{para: "-t must > 0".to_string()})
                }
                thread
            },
            None => 4_usize,
        },
        outpath: match para.outpath { // 输出结果路径，不存在则自动创建，默认./
            Some(o) => {
                let tmp_path = PathBuf::from(&o);
                if !(tmp_path.exists() && tmp_path.is_dir()) {
                    if let Err(err) = create_dir_all(&tmp_path) {
                        return Err(MyError::CreateDirAllError{dir_name: o, error: err})
                    }
                }
                tmp_path
            },
            None => PathBuf::from("./"),
        },
    };
    Ok(out)
}
