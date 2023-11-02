/*
 * Copyright (c) 2022, Beatrice Borsari
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


// Define parameters

params.maf = 0.01
params.r2 = 0.30
params.l = 200000
params.index = null
params.outFolder = 'results'
params.help = false


// Print usage

if (params.help) {
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run filter.imputed.variants.nf [options]'
    log.info ''
    log.info 'Parameters:'
    log.info ' --maf MAF			MAF threshold (default = 0.01)'
    log.info ' --r2  R-SQUARED  	        R2 threshold (default = 0.3)'
    log.info ' --l VARIANTS/CHUNK       	# of variants per chunk (default: 200000)'
    log.info ' --index INDEX             	index file containing name of chr and corresponding vcf file'
    log.info ' --outFolder DIRECTORY		output directory (default: results)'
    log.info ''
    exit(1)
}


// Check mandatory parameters

if (!params.index) {
    exit 1, "Index file not specified."
}


// Print parameter selection

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "MAF threshold               : ${params.maf}"
log.info "R2 threshold                : ${params.r2}"
log.info "Variants/chunk              : ${params.l}"
log.info "Index file                  : ${params.index}"
log.info "Output directory            : ${params.outFolder}"
log.info ''


//~~~~~~~~~~~~~~
// BEGIN 
//~~~~~~~~~~~~~~


// init channels for chrom and vcf files
index = file(params.index)

// create two identical channels that contain
// 1. name of chr
// 2. vcf of chr
Channel.from(index.readLines())
.map { line ->
  def list = line.tokenize()
  def chr = list[0]
  def vcf = resolveFile(list[1], index)
  [ chr, vcf ]
}.into
{chr_vcf_1_ch;chr_vcf_2_ch;chr_vcf_4_ch} 


// add, to one of the channels above,
// the name of chr index obtained with tabix
chr_vcf_2_ch.map{it -> [it[0], it[1], it[1]+".tbi"]}.set{chr_vcf_3_ch}


// for each chr, split VCF 
// in chuncks each containing n positions
process split_vcf {

  input:
  tuple chr, file(vcf) from chr_vcf_1_ch

  output:
  tuple chr, file("chunk*") into chunks_ch    

  script:
  """
  bcftools query -f '%CHROM\t%POS\n' $vcf > positions
  split -d -a 10 -l ${params.l} positions chunk
  """
}

// create a channel that contains
// 1. chr
// 2. vcf
// 3. index
// 4. chunk 
// 5. number of chunk
chr_vcf_3_ch.combine(chunks_ch.transpose(), by: 0)
.map{it -> [it[0], it[1], it[2], it[3], it[3].baseName]}
.set{ready2parse_chunks_ch}

// filter vcf file based on MAF and R2 
process filter_vcf {

  input:
  tuple chr, file(vcf), file(index), file(chunk), n_chunk from ready2parse_chunks_ch

  output:
  tuple chr, file("${n_chunk}.filtered.vcf") into filtered_chunks_ch

  script:
  """
  start=\$(head -1 ${chunk} | cut -f2)
  end=\$(tail -1 ${chunk} | cut -f2)

  tabix $vcf ${chr}":"\${start}"-"\${end} > sub.vcf
  
  parse.py -g sub.vcf -m ${params.maf} -r ${params.r2} -o ${n_chunk}.filtered.vcf
  """
}


// create a channel that contains:
// 1. chr
// 2. original vcf (needed to extract the header)
// 3. all filtered chunks of vcf files
// chunks are sorted by chunk number
chr_vcf_4_ch.combine(filtered_chunks_ch.groupTuple(by : 0, sort: {it.name}), by: 0)
.set{ready2merge_chunks_ch}


// merge all filtered chunks of a chromosome
// in one file
process merge_filtered_chunks {

  publishDir "${params.outFolder}/${chr}/filtered", overwrite: true

  input:
  tuple chr, vcf, file(filtered_chunks) from ready2merge_chunks_ch

  output:
  tuple chr, file("${chr}.filtered.vcf.gz"), file("${chr}.filtered.infoColumn.tsv") into chr_filtered_vcf_ch

  script:
  """
  # extract header from original vcf
  zcat $vcf | head -19 > "${chr}.filtered.vcf"
 
  # create header for infoColumn tsv file
  echo -e "snp_id\tAF\tMAF\tR2\tER2\tsnp_type" > ${chr}.filtered.infoColumn.tsv

  # merge filtered chunks
  # and extract informative columns
  for f in ${filtered_chunks}; do
        if [ -s \$f ];
        then
	    ## save content of chunk
	    cat \$f >> "${chr}.filtered.vcf"

            ## extract snp and info columns
            cut -f3,8 \$f | awk 'BEGIN{FS=OFS="\t"}
				      {n=split(\$2, a, ";");
                                       for (i=1;i<n;i++) {split(a[i], b, "="); a[i]=b[2]}; 
                                       if (n<5){print \$1, a[1], a[2], a[3], "NA", a[4]} else {print \$1, a[1], a[2], a[3], a[4], a[5]
                                      }}' >> "${chr}.filtered.infoColumn.tsv"

        fi 
  done

  # compress file
  bgzip "${chr}.filtered.vcf"
  """

}
 
chr_filtered_vcf_ch.view()



/*
 * Given a string path resolve it against the index file location.
 * Params:
 * - str: a string value represting the file path to be resolved
 * - index: path location against which relative paths need to be resolved
 */
def resolveFile( str, index ) {
  if( str.startsWith('/') || str =~ /^[\w\d]*:\// ) {
    return file(str)
  }
  else if( index instanceof Path ) {
    return index.parent.resolve(str)
  }
  else {
    return file(str)
  }
}

