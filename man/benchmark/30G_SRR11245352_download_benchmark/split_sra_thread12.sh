#!/bin/bash
#SBATCH -J "split12"
#SBATCH -p Acluster
#SBATCH --exclusive
#SBATCH -n 12
#SBATCH --output=%j.out
#SBATCH --error=%j.err

export PATH=/BioII/wangjb/20T/songyabing/software/sratoolkit.3.1.1-centos_linux64/bin:/BioII/wangjb/20T/songyabing/software/parallel-fastq-dump-master:$PATH

# create time folder
time_folder=/BioII/wangjb/20T/songyabing/benchfastq/SRR11245352_time
# mkdir -p ${time_folder}
for tr in 12
  do
  echo -e "tool\ttype\ttr\ttime\tbatch" >${time_folder}/SRA_time_${tr}.txt
  echo -e "tool\ttype\ttr\ttime\tbatch" >${time_folder}/ENA_time_${tr}.txt
  for ((i=1; i<=3; i++))
    do
    # ======================================= split sra into fastq =============================================
    # prepare folder
    mkdir -p /BioII/wangjb/20T/songyabing/benchfastq/tr_${tr}
    cp -r /BioII/wangjb/20T/songyabing/benchfastq/SRA_2/* /BioII/wangjb/20T/songyabing/benchfastq/tr_${tr}
    cd /BioII/wangjb/20T/songyabing/benchfastq/tr_${tr}

    if [ $tr -eq 1 ]
    then
      # fastq-dump does not support multithreading
      # split sra to fastq with fastq-dump (gzip)
      start=$(date +%s)
      fastq-dump --split-files -O ./fastq_dump_gz --gzip ./SRR11245352/SRR11245352.sra
      end=$(date +%s)
      fastq_dump_gz_time=$((end - start))
      echo -e "fastq-dump\tgz\t${tr}\t${fastq_dump_gz_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt
      rm -rf ./fastq_dump_gz

      # split sra to fastq with fastq-dump (ungzip)
      start=$(date +%s)
      fastq-dump --split-files -O ./fastq_dump_fq ./SRR11245352/SRR11245352.sra
      end=$(date +%s)
      fastq_dump_fq_time=$((end - start))
      echo -e "fastq-dump\tfq\t${tr}\t${fastq_dump_fq_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt
      rm -rf ./fastq_dump_fq
    fi

    # fasterq-dump will ignore --threads 1: https://github.com/ncbi/sra-tools/issues/494
    if [ $tr -gt 1 ]
    then
      # split sra to fastq with fasterq-dump (pigz)
      start=$(date +%s)
      fasterq-dump --split-files -O ./fasterq_dump_pigz --threads ${tr} ./SRR11245352/SRR11245352.sra
      pigz -p ${tr} ./fasterq_dump_pigz/*fastq
      end=$(date +%s)
      fasterq_dump_pigz_time=$((end - start))
      echo -e "fasterq-dump\tpigz\t${tr}\t${fasterq_dump_pigz_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt
      rm -rf ./fasterq_dump_pigz

      # split sra to fastq with fasterq-dump (gzip)
      start=$(date +%s)
      fasterq-dump --split-files -O ./fasterq_dump_gzip --threads ${tr} ./SRR11245352/SRR11245352.sra
      gzip ./fasterq_dump_gzip/*fastq
      end=$(date +%s)
      fasterq_dump_gzip_time=$((end - start))
      echo -e "fasterq-dump\tgzip\t${tr}\t${fasterq_dump_gzip_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt
      rm -rf ./fasterq_dump_gzip

      # split sra to fastq with fasterq-dump
      start=$(date +%s)
      fasterq-dump --split-files -O ./fasterq_dump_fq --threads ${tr} ./SRR11245352/SRR11245352.sra
      end=$(date +%s)
      fasterq_dump_fq_time=$((end - start))
      echo -e "fasterq-dump\tfq\t${tr}\t${fasterq_dump_fq_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt
      rm -rf ./fasterq_dump_fq
    fi

    # split sra to fastq with parallel-fastq-dump (gz)
    start=$(date +%s)
    parallel-fastq-dump --split-files -O ./parallel_fastq_dump_gz --threads ${tr} --gzip --sra-id ./SRR11245352/SRR11245352.sra
    end=$(date +%s)
    parallel_fastq_dump_gz_time=$((end - start))
    echo -e "parallel-fastq-dump\tgz\t${tr}\t${parallel_fastq_dump_gz_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt
    rm -rf ./parallel_fastq_dump_gz

    # split sra to fastq with parallel-fastq-dump (fq)
    start=$(date +%s)
    parallel-fastq-dump --split-files -O ./parallel_fastq_dump_fq --threads ${tr} --tmpdir ./ --sra-id ./SRR11245352/SRR11245352.sra
    end=$(date +%s)
    parallel_fastq_dump_fq_time=$((end - start))
    echo -e "parallel-fastq-dump\tfq\t${tr}\t${parallel_fastq_dump_fq_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt
    rm -rf ./parallel_fastq_dump_fq
  done
  rm -rf /BioII/wangjb/20T/songyabing/benchfastq/tr_${tr}
done

