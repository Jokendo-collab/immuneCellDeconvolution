
setwd("C:/Users/Javan_Okendo/Desktop/cybersort/cybersort_source_code")

dir()
source("CIBERSORT.R")

#Run the individual signature matrix file from each patient groups
ltbi_result <- CIBERSORT("LTBI_signature_matrix","LM22.txt",perm = 100,QN = TRUE,absolute = FALSE,abs_method = 'sig.score')

recTB_result <- CIBERSORT("reccurentTB_signature_matrix","LM22.txt",perm = 1000,QN = TRUE,absolute = FALSE,abs_method = 'sig.score')
          

prevTB_result <- CIBERSORT("prevTB_signature_matrix","LM22.txt",perm = 1000, QN = TRUE,absolute = FALSE,abs_method = 'sig.score')
