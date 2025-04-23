use std::path::PathBuf;
use std::time::SystemTime;
use std::thread;

/// parse_para: 解析命令行参数
/// MyError: 自定义错误类型
use mutate::{
    parse_paras::parse_para,
    error::MyError,
    mutation::mutate_fastq,
    utils::{split_2_group, my_time},
};

fn main() {
    if let Err(e) = run() {
        println!("{}", e); // 这里不要用`{:?}`，会打印结构体而不是打印指定的错误信息
    }
}

fn run() -> Result<(), MyError> {
    // 开始时间
    let t1 = SystemTime::now();

    // 解析参数
    let paras = parse_para()?;

    // 将所有要突变的fastq文件分为thread部分，使用多线程进行突变
    thread::scope(|s| -> Result<(), MyError> {
        // 使用scope调用线程
        // 解决多个常规thread::spawn读取数组的生命周期问题：
        //     https://stackoverflow.com/questions/75000029/reading-a-vector-from-multiple-threads
        //     https://stackoverflow.com/questions/32750829/how-can-i-pass-a-reference-to-a-stack-variable-to-a-thread
        //     https://stackoverflow.com/questions/73950960/proper-way-to-share-references-to-vec-between-threads
        // scope与channel联合使用：https://blog.logrocket.com/using-rust-scoped-threads-improve-efficiency-safety/
        // scope内for循环调用线程，进队指定变量进行move：
        //     https://stackoverflow.com/questions/73747874/why-is-threadscope-unable-to-copy-primitive-types-without-a-move-closure
        //     https://stackoverflow.com/questions/75126641/value-moved-into-closure-in-previous-iteration-of-loop
        let group = split_2_group(paras.thread, paras.fastq.len());
        for g in group.iter() { // 调用指定数量线程
            let start = g[0];
            let end = g[1];
            s.spawn({ // 这里在闭包外又包裹了一层{}，以实现move时移动的是all_fastq的引用
                let tmp_fastq = &paras.fastq[start..end];
                let tmp_outpath = paras.outpath.clone();
                move || -> Result<(), MyError> {
                    mutate_group_fastq(tmp_fastq, tmp_outpath, paras.ratio, paras.pos)
                }
            });
        }
        Ok(())
    })?;

    // 结束时间
    let t2 = SystemTime::now();

    // 打印总耗时
    println!("{}", my_time("[Total time]", t2.duration_since(t1).unwrap().as_nanos()));

    Ok(())
}

/// 遍历每个fastq文件进行突变
fn mutate_group_fastq(fastq_files: &[PathBuf], outpath: PathBuf, ratio: f32, save_pos: bool) -> Result<(), MyError> {
    for f in fastq_files {
        mutate_fastq(&f, &outpath, ratio, save_pos)?;
    }
    Ok(())
}
