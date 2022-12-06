
x <- readLines('u5_100.hmm')
write(c('>u5', paste0(toupper(sub('\\s\\-', '', unlist(stringr::str_extract_all(x, '[atcg]\\s\\-')))), collapse = '')), 
      file = 'u5_100.model.fasta')

x <- readLines('u3_100_rc.hmm')
write(c('>u3', paste0(toupper(sub('\\s\\-', '', unlist(stringr::str_extract_all(x, '[atcg]\\s\\-')))), collapse = '')), 
      file = 'u3_100.model.fasta')