#!/bin/bash

threads=10
opt_minimap2=""


set -ue
tmpdir=$(mktemp -d -p ${TMPDIR:-/tmp})
trap "[[ $tmpdir ]] && rm -rf $tmpdir" 0 1 2 3 15



function usage()
{
  cat <<EOF

Genome-guided error correction for error prone sequencing of cDNA amplicon



Usage
---

To perform all steps of error correction, please use 'correct' subcommand as:

  $0 correct -i INFILE.fq  -g GENOME [-p PASS(default:2, or 1)] [-t THREADS(default:10)] [options used in individual steps]


Each of the individual steps can be run by other subcommands,
'ggc' or 'ecc', as:

  # genome guided clustering (via BLAT)
  $0 ggc -i INFILE -g GENOME [-t THREADS(default:10)] [-l SEQLEN_MIN(default:1000, should not be changed for Canu error correction)]  [-e error_correction_seqnum_min(default:5)]

  # error correction within each of the clusters (with Canu)
  $0 ecc -i INFILE -c CLUSTER_FILE('ggc' output) [-t THREADS(default:10)]



Note that sequence IDs should not include comma (',').


How it works
------------
Reads are grouped by isoform based on their alignments with
the specified genome sequence (by minimap2), and the orphan
reads (that does not grouped with others) are included
to the grouop that has the maximum length of exons in total.
Sequence error correction is performed within each of the
groups.

For two pass mode, the corrected sequences (in the first pass)
are used for grouping, and sequence error correction itself is
performed on the uncorrected reads.

Note that Canu does not perform error correction on the sequences shorter
than 1000bp (in default), and they are discarded at the beginning.



Requirements
------------
Software below is required to run "$0". Their versions used for
development are indicated in '( )'

* minimap2 (version 2.12)
* setqk (version 1.3)
* bedtools [bed12ToBed6, mergeBed, groupBy] (version 2.27.1)
* samtools (version 1.7)
* Canu (version 1.7.1)
* R (with package 'stringr')
* GNU awk (ver.4.1.3)



Author
------
This software is written by Hideya Kawaji.
Copyright (c) 2018-2020 RIKEN & Tokyo Metropolitan Institute of Medical Science.
Distributed under Apache License 2.0.


EOF
  exit 1;
}



function _incorporate_tiny_groups()
{
  local infile=$1
  local error_correction_seqnum_min=$2

  ### incorporate tiny groups into regular ones.
  cat <<EOF | R --slave
  library(stringr)
  length_common_exons <- function( a, b )
  {
    ce = intersect( unlist( str_split(a,",") ) , unlist(str_split(a,",")) )
    if ( length(ce) == 0 ){ return(0) }
    coord = str_split(ce,"_", n=4, simplify=T)
    sum( as.integer( coord[,3] ) - as.integer( coord[,2] ) )
  }
  max_common_exons <- function( a, bVec )
  {
    if ( length( bVec) == 0 ) { return(NA) }
    common_length = sapply(bVec, function(b) length_common_exons(a, b))
    if ( max(common_length) == 0 ) { return(NA) }
    bVec[ which.max( common_length ) ]
  }
  
  ggc = read.table("${infile}", sep="\t", as.is=T, row.names=1)
  n_reads = str_count(unlist(ggc),",")
  idx = which( n_reads >= ${error_correction_seqnum_min} )
  if ( (length(idx) == 0) | (length(idx) == nrow(ggc)) ) { return() }
  cluster = rownames(ggc)[idx]
  orphan = rownames(ggc)[-idx]
  for (o in orphan)
  { 
    clst = max_common_exons(o, cluster)
    if ( ! is.na(clst) )
    { 
      ggc[clst,1] = sprintf("%s,%s", ggc[clst,1], ggc[o,1])
    }
  }
  write.table( cbind( "#reads" = ggc[cluster,1], "exons" = cluster ) , sep="\t", quote=F, row.names=F, col.names=F )
EOF
}



function _error()
{
  printf "Error: %s\n" $1
  exit 1
}



function align()
{
  local infile=
  local genome=
  local seqlen_min=1000 # should not be changed for Canu's correction (since Canu won't correct shorter sequences in default)
  local threads=${threads}
  local OPTIND_OLD=$OPTIND
  OPTIND=1
  while getopts i:g:t:l: opt
  do
    case ${opt} in
    i) infile=${OPTARG};;
    g) genome=${OPTARG};;
    t) threads=${OPTARG};;
    l) seqlen_min=${OPTARG};;
    *) _error "align" ;;
    esac
  done
  OPTIND=${OPTIND_OLD}
  if [ ! -n "${infile-}" ]; then _error "align"; fi
  if [ ! -n "${genome-}" ]; then _error "align"; fi
  if [ ! -n "${seqlen_min-}" ]; then _error "align"; fi

  #local infile=$1
  #local genome=$2
  #local seqlen_min=$3
  local tmpdir=$(mktemp -d -p ${tmpdir})

  ### prep
  seqtk comp ${infile} \
  | awk --assign seqlen_min=$seqlen_min '{if ($2 >= seqlen_min){ print $1} }' \
  > ${tmpdir}/infile.seqids

  seqtk subseq ${infile} ${tmpdir}/infile.seqids \
  > ${tmpdir}/infile

  ### align
  minimap2 -t ${threads} -ax splice ${genome} ${opt_minimap2} ${tmpdir}/infile \
  > ${tmpdir}/infile.sam

  samtools view -uSF 0x800 ${tmpdir}/infile.sam  \
  | bamToBed -bed12 \
  | sort -k1,1 -k2,2n
}



function ggc()
{
  local error_correction_seqnum_min=$1

  # take bed12 at stdin
  local tmpdir=$(mktemp -d -p ${tmpdir})

  ### group by exon_patterns
  bed12ToBed6 \
  | sort -k1,1 -k2,2n \
  > ${tmpdir}/infile.bed6

  cat ${tmpdir}/infile.bed6 \
  | mergeBed -i stdin \
  | awk '{printf "%s\t%s\t%s\t%s_%s_%s\n", $1,$2,$3,$1,$2,$3}' \
  | intersectBed -wa -wb -a - -b ${tmpdir}/infile.bed6  \
  | awk '{printf "%s\t%s\n", $8,$4}' \
  | uniq | sort -k1,1 | groupBy -i - -grp 1 -o distinct -opCols 2 \
  | sort -k2,2 | groupBy -i - -grp 2 -o distinct -opCols 1 \
  > ${tmpdir}/infile.cluster

  _incorporate_tiny_groups ${tmpdir}/infile.cluster $error_correction_seqnum_min
}



function ecc()
{
  local infile=
  local infile_cluster=

  local OPTIND_OLD=$OPTIND
  OPTIND=1
  while getopts i:t:c: opt
  do
    case ${opt} in
    i) infile=${OPTARG};;
    c) infile_cluster=${OPTARG};;
    t) threads=${OPTARG};;
    *) echo "Error in ecc\n" 1>&2; exit 1;;
    esac
  done
  OPTIND=${OPTIND_OLD}
  if [ ! -n "${infile-}" ]; then _error "ecc"; fi
  if [ ! -n "${infile_cluster-}" ]; then _error "ecc"; fi
  if [ ! -n "${threads-}" ]; then _error "ecc"; fi
  local tmpdir=$(mktemp -d -p ${tmpdir})

  ### prep
  cat ${infile_cluster} \
  | awk --assign prefix=${tmpdir} '{
      if ( match($0,"^#") ){ next }
      outfile=sprintf("%s/%010X.seqids",prefix,NR);
      len=gsub(",","\n",$1);
      print $1 > outfile;
      close(outfile)
    }'

  ### correction for each cluster
  find ${tmpdir} -name '*.seqids' \
  | xargs --max-procs=${threads} -L 1 -I {} bash -c \
    "mkdir -p {}.canu ;
     seqtk subseq ${infile} {} | gzip -c > {}.canu/clst.fq.gz ; 
     cd {}.canu ;
     len_k=\$(( \$(seqtk comp clst.fq.gz | cut -f 2 |sort -nr |head -1) / 1000 + 1 ))k;
     canu -nanopore-raw clst.fq.gz -correct -p clst.canu genomeSize=\${len_k} corOutCoverage=999 corMinCoverage=0 stopOnReadQuality=false > /dev/null ;
     pwd 1>&2
     "
  find ${tmpdir} -name '*.correctedReads.fasta.gz' \
  | xargs -L 1 -I {} bash -c \
    "gunzip -c {} | awk --assign bn=\$(basename \$(dirname {} .canu)) '{ if(match(\$0,/^>/)){print \$0,\"correction_group:\"bn}else{print \$0}  }'"
}



function correct()
{ 
  local infile=
  local genome=
  local pass=2
  local seqlen_min=1000 # should not be changed for Canu's correction (since Canu won't touch shorter sequences)
  local error_correction_seqnum_min=5
  local outfile_onepass=
  local outfile_twopass=

  while getopts i:g:t:l:p:1:2: opt
  do
    case ${opt} in
    i) infile=${OPTARG};;
    g) genome=${OPTARG};;
    t) threads=${OPTARG};;
    l) seqlen_min=${OPTARG};;
    p) pass=${OPTARG};;
    1) outfile_onepass=${OPTARG};;
    2) outfile_twopass=${OPTARG};;
    *) _error "correct";;
    esac
  done
  if [ ! -n "${pass-}" ]; then _error "correct"; fi
  if [ ! -n "${infile-}" ]; then _error "correct"; fi
  if [ ! -n "${genome-}" ]; then _error "correct"; fi
  if [ ! -n "${threads-}" ]; then _error "correct"; fi
  local tmpdir=$(mktemp -d -p ${tmpdir})

  ### align
  align -i ${infile} -g ${genome} -t $threads -l $seqlen_min \
  > ${tmpdir}/infile.bed12

  ### clustering
  cat ${tmpdir}/infile.bed12 \
  | ggc $error_correction_seqnum_min \
  > ${tmpdir}/ggc.out

  ### error correction
  ecc -i ${infile} -c ${tmpdir}/ggc.out -t $threads \
  | tee ${tmpdir}/ecc.out > ${tmpdir}/onepass.out

  if [ -n "${outfile_onepass-}" ]; then
    cp -f ${tmpdir}/ggc.out $outfile_onepass
  fi

  case $pass in
  1)
    cat ${tmpdir}/onepass.out
    ;;
  2)
    ### realign
    align -i ${tmpdir}/onepass.out -g ${genome} -t $threads -l $seqlen_min \
    > ${tmpdir}/onepass.bed12

    ### obtain the initial alignments for the uncorrected reads
    cut -f 4 ${tmpdir}/onepass.bed12 \
    | sort > ${tmpdir}/onepass.bed12.id

    join -v 1 -t "	" -1 4 -2 1 \
      <( sort -k4,4 ${tmpdir}/infile.bed12 ) \
      ${tmpdir}/onepass.bed12.id \
    | awk 'BEGIN{OFS="\t"}{n=$1;$1=$2;$2=$3;$3=$4;$4=n;print}' \
    | sort -k1,1 -k2,2n \
    >> ${tmpdir}/onepass.bed12

    ### clustering
    cat ${tmpdir}/onepass.bed12 \
    | ggc $error_correction_seqnum_min \
    > ${tmpdir}/onepass.ggc.out

    ### error correction
    ecc -i ${infile} -c ${tmpdir}/onepass.ggc.out -t $threads \
    | tee ${tmpdir}/onepass.ecc.out > ${tmpdir}/twopass.out

    if [ -n "${outfile_twopass-}" ]; then
      cp -f ${tmpdir}/onepass.ggc.out $outfile_twopass
    fi

    cat ${tmpdir}/twopass.out
    ;;
  *)
    _error "correct"
    ;;
  esac
}


###      ###
### main ###
###      ###
if [ $# == 0 ] ; then usage; fi
subcommand="$1"
shift
case $subcommand in
  ggc) ggc $@;;
  ecc) ecc $@;;
  correct) correct $@;;
  *) usage;;
esac


