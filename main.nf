params.rawdata = null
params.publishDir = "/tmp"
params.dockerImage = 'geneplaza/sapda-k5:debian2'

rawFileChannel = Channel.fromPath(params.rawdata)

process inputFile {

  echo true
  container params.dockerImage

  input:
  file rawdata from rawFileChannel
  output:
  file("TEST.txt") into inputFileOutChan

  script:
    """
    unzip *.zip
    ls -l
    for i in {1..22}; do
      cat id_*_chr\${i}.23andme.txt >> result.txt;
    done
    awk -F '\t' 'NR==FNR{c[\$1]++;next};c[\$1]' OFS="\t" /reference-file/Ancestry_59K_SNPs.txt result.txt > TEST.txt
    #if [ -f id_*_chrX.23andme.txt ]; then
    #  cat id_*_chrX.23andme.txt >> result.txt;
    #fi
    """
}

process app {

  publishDir params.publishDir, mode: 'copy'
  echo true
  container params.dockerImage

  input:
  file rawdata from inputFileOutChan

  output:
  file ("2nd_model/DEEP/gnuplot/*.dat") into appOut
  file ("ADMIX_OLD_RECENT_MEAN6.txt") into appOut2
  file ("GSI_OLD_RECENT_MEAN6.txt") into appOut3
  file ("*") into appOut4

  script:
  """
  echo \$PWD
  cp ${rawdata} /app/TEST.txt
  pushd /app
  /app/SCRIPT_5POP_ALT_GenePlaza.sh
  popd
  mv /app/* .
  """
}

process transform {
  publishDir params.publishDir, mode: 'copy'
  echo true
  container 'amancevice/pandas:slim'

  input:
  file('*') from appOut.collect()
  file admixFile from appOut2
  file gsiFile from appOut3

  output:
  file 'result.json' into transformOut

  """
  #!/usr/local/bin/python

  import pandas as pd
  import os
  import glob
  import ast
  import json

  json_out = {}

  admixture_path = 'ADMIX_OLD_RECENT_MEAN6.txt'
  df_admix = pd.read_csv(admixture_path, delimiter='\t')
  json_out['admixture-proportions'] = ast.literal_eval(df_admix.to_json(orient='records'))

  gsi_path = 'GSI_OLD_RECENT_MEAN6.txt'
  df_gsi = pd.read_csv(gsi_path, delimiter='\t')
  json_out['gsi-proportions'] = ast.literal_eval(df_gsi.to_json(orient='records'))
  shared_genetic_drift_paths = glob.glob('*.dat')

  for path in shared_genetic_drift_paths:
      population = (os.path.basename(path.strip(".dat")))
      df_k = pd.read_csv(path, delimiter='\t', header=None)
      df_k = df_k.drop(df_k.columns[range(13,24)], axis=1)
      df_k['Chromosome'], df_k['RSID'] = df_k[0].str.split('-', 1).str
      df_k = df_k.drop(df_k.columns[0], axis=1)
      df_k['genotype'] = df_k[[1, 2]].apply(lambda x: ''.join(x), axis=1)
      cols = [0,1]
      df_k = df_k.drop(df_k.columns[cols], axis =1)
      df_k.columns=['reference', 'E-Asian-AF','SE-Asian-AF','Siberian-Amerindian-AF','S-Asian-AF','NE-European-AF','W-European-AF','S-European-AF','W-African-AF','E-African-AF','chromosome','RSID','genotype']
      json_out[population] = ast.literal_eval(df_k.to_json(orient='records'))

  with open('result.json', 'w') as outfile:
    json.dump(json_out, outfile, ensure_ascii=False)
  """
}


 workflow.onError {
   errFile = file("${params.publishDir}/errorMsg.txt")
   errFile.write("${workflow.errorMessage}")
}
