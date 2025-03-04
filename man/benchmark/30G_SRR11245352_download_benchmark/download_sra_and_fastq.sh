#!/bin/bash

# create time folder
time_folder=/home/syb/projects/04_scfetch/benchfastq/SRR11245352_time
mkdir -p ${time_folder}
# output file
echo -e "tool\ttype\ttime\tbatch" >${time_folder}/SRA_time_down.txt
echo -e "tool\ttype\ttime\tbatch" >${time_folder}/ENA_time_down.txt
for ((i=1; i<=10; i++))
  do
  # ======================================= download sra from SRA =============================================
  # prepare folder
  cd /home/syb/projects/04_scfetch/benchfastq
  mkdir SRA_${i}
  cd SRA_${i}
  # download sra from SRA
  start=$(date +%s)
  prefetch -X 100G SRR11245352
  end=$(date +%s)
  prefetch_time=$((end - start))
  echo -e "prefetch\tsra\t${tr}\t${prefetch_time}\t${i}" >>${time_folder}/SRA_time_${tr}.txt

  # ======================================= download sra from ENA =============================================
  # prepare folder
  cd /home/syb/projects/04_scfetch/benchfastq
  mkdir ENA_${i}
  cd ENA_${i}
  # download sra from ENA
  start=$(date +%s)
  ascp -QT -l 300m -P33001 -k 1 -i /home/syb/miniconda3/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/srr/SRR112/052/SRR11245352 .
  mv SRR11245352 SRR11245352.sra
  end=$(date +%s)
  ascp_time=$((end - start))
  echo -e "ascp\tsra\t${tr}\t${ascp_time}\t${i}" >>${time_folder}/ENA_time_${tr}.txt

  # ======================================= download fastq from ENA =============================================
  # prepare folder
  cd /home/syb/projects/04_scfetch/benchfastq/ENA_${i}
  # download sra from ENA
  start=$(date +%s)
  ascp -QT -l 300m -P33001 -k 1 -i /home/syb/miniconda3/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR112/052/SRR11245352/SRR11245352_1.fastq.gz .
  ascp -QT -l 300m -P33001 -k 1 -i /home/syb/miniconda3/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR112/052/SRR11245352/SRR11245352_2.fastq.gz .
  end=$(date +%s)
  ascp_fastq_time=$((end - start))
  echo -e "ascp_fastq\tfq\t${tr}\t${ascp_fastq_time}\t${i}" >>${time_folder}/ENA_time_${tr}.txt
  rm -rf /home/syb/projects/04_scfetch/benchfastq/ENA_${i}
done
