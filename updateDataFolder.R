AAVengeR_dir <- '~/AAVengeR'

system(paste0('wget -q -O ', AAVengeR_dir, '/data/blastDBs.tar https://microb120.med.upenn.edu/data/AAVengeR/blatDBs.tar'))
system(paste0('tar xf ', AAVengeR_dir, '/data/blastDBs.tar -C ', AAVengeR_dir, '/data/'))
system(paste0('rm -f ', AAVengeR_dir, '/data/blastDBs.tar'))

system(paste0('wget -q -O ', AAVengeR_dir, '/data/genomeAnnotations.tar https://microb120.med.upenn.edu/data/AAVengeR/genomeAnnotations.tar'))
system(paste0('tar xf ', AAVengeR_dir, '/data/genomeAnnotations.tar -C ', AAVengeR_dir, '/data/'))
system(paste0('rm -f ', AAVengeR_dir, '/data/genomeAnnotations.tar'))