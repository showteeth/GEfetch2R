library(tidyverse)
library(GEfetch2R)

# Count matrix -------------
# * Smart-Seq2 ---------------
GSE94820.time = system.time({GSE94820.cnt = ParseGEO(acce = "GSE94820", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE94820.cnt[1:5,1:5]
GSE214129.time = system.time({GSE214129.cnt = ParseGEO(acce = "GSE214129", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE214129.cnt[1:5,1:5]
GSE292261.time = system.time({GSE292261.cnt = ParseGEO(acce = "GSE292261", down.supp = TRUE, supp.idx = 2, supp.type = "count", load2R = F)})
# GSE292261.cnt[1:5,1:5]
GSE224915.time = system.time({GSE224915.cnt = ParseGEO(acce = "GSE224915", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE224915.cnt[1:5,1:5]
GSE185076.time = system.time({GSE185076.cnt = ParseGEO(acce = "GSE185076", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE185076.cnt[1:5,1:5]
GSE241004.time = system.time({GSE241004.cnt = ParseGEO(acce = "GSE241004", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE241004.cnt[1:5,1:5]
GSE244538.time = system.time({GSE244538.cnt = ParseGEO(acce = "GSE244538", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE244538.cnt[1:5,1:5]
GSE266889.time = system.time({GSE266889.cnt = ParseGEO(acce = "GSE266889", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE266889.cnt[1:5,1:5]
GSE297431.time = system.time({GSE297431.cnt = ParseGEO(acce = "GSE297431", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE297431.cnt[1:5,1:5]
GSE286046.time = system.time({GSE286046.cnt = ParseGEO(acce = "GSE286046", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE286046.cnt[1:5,1:5]
GSE262789.time = system.time({GSE262789.cnt = ParseGEO(acce = "GSE262789", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE262789.cnt[1:5,1:5]
GSE276332.time = system.time({GSE276332.cnt = ParseGEO(acce = "GSE276332", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE276332.cnt[1:5,1:5]
GSE274944.time = system.time({GSE274944.cnt = ParseGEO(acce = "GSE274944", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE274944.cnt[1:5,1:5]
GSE266915.time = system.time({GSE266915.cnt = ParseGEO(acce = "GSE266915", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE266915.cnt[1:5,1:5]
GSE182219.time = system.time({GSE182219.cnt = ParseGEO(acce = "GSE182219", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE182219.cnt[1:5,1:5]
GSE225056.time = system.time({GSE225056.cnt = ParseGEO(acce = "GSE225056", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE225056.cnt[1:5,1:5]
GSE280276.time = system.time({GSE280276.cnt = ParseGEO(acce = "GSE280276", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE280276.cnt[1:5,1:5]
GSE268201.time = system.time({GSE268201.cnt = ParseGEO(acce = "GSE268201", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)}) # with comment char: '#
# GSE268201.cnt[1:5,1:5]
GSE252076.time = system.time({GSE252076.cnt = ParseGEO(acce = "GSE252076", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE252076.cnt[1:5,1:5]
GSE284380.time = system.time({GSE284380.cnt = ParseGEO(acce = "GSE284380", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE284380.cnt[1:5,1:5]
GSE245172.time = system.time({GSE245172.cnt = ParseGEO(acce = "GSE245172", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE245172.cnt[1:5,1:5]
GSE215960.time = system.time({GSE215960.cnt = ParseGEO(acce = "GSE215960", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE215960.cnt[1:5,1:5]
GSE261729.time = system.time({GSE261729.cnt = ParseGEO(acce = "GSE261729", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE261729.cnt[1:5,1:5]
GSE266449.time = system.time({GSE266449.cnt = ParseGEO(acce = "GSE266449", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE266449.cnt[1:5,1:5]
GSE264105.time = system.time({GSE264105.cnt = ParseGEO(acce = "GSE264105", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)}) # anomalous data
# GSE264105.cnt[1:5,1:5]
GSE249746.time = system.time({GSE249746.cnt = ParseGEO(acce = "GSE249746", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE249746.cnt[1:5,1:5]
GSE268974.time = system.time({GSE268974.cnt = ParseGEO(acce = "GSE268974", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE268974.cnt[1:5,1:5]
GSE241165.time = system.time({GSE241165.cnt = ParseGEO(acce = "GSE241165", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE241165.cnt[1:5,1:5]
GSE271100.time = system.time({GSE271100.cnt = ParseGEO(acce = "GSE271100", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE271100.cnt[1:5,1:5]
GSE264553.time = system.time({GSE264553.cnt = ParseGEO(acce = "GSE264553", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE264553.cnt[1:5,1:5]
GSE267774.time = system.time({GSE267774.cnt = ParseGEO(acce = "GSE267774", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE267774.cnt[1:5,1:5]
GSE261258.time = system.time({GSE261258.cnt = ParseGEO(acce = "GSE261258", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE261258.cnt[1:5,1:5]
GSE268038.time = system.time({GSE268038.cnt = ParseGEO(acce = "GSE268038", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE268038.cnt[1:5,1:5]
GSE234572.time = system.time({GSE234572.cnt = ParseGEO(acce = "GSE234572", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE234572.cnt[1:5,1:5]
GSE265909.time = system.time({GSE265909.cnt = ParseGEO(acce = "GSE265909", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE265909.cnt[1:5,1:5]
GSE237867.time = system.time({GSE237867.cnt = ParseGEO(acce = "GSE237867", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE237867.cnt[1:5,1:5]
GSE256466.time = system.time({GSE256466.cnt = ParseGEO(acce = "GSE256466", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE256466.cnt[1:5,1:5]
GSE211550.time = system.time({GSE211550.cnt = ParseGEO(acce = "GSE211550", down.supp = TRUE, supp.idx = 2, supp.type = "count", load2R = F)})
# GSE211550.cnt[1:5,1:5]
GSE207479.time = system.time({GSE207479.cnt = ParseGEO(acce = "GSE207479", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE207479.cnt[1:5,1:5]
GSE241089.time = system.time({GSE241089.cnt = ParseGEO(acce = "GSE241089", down.supp = TRUE, supp.idx = 5, supp.type = "count", load2R = F)})
# GSE241089.cnt[1:5,1:5]
GSE230792.time = system.time({GSE230792.cnt = ParseGEO(acce = "GSE230792", down.supp = TRUE, supp.idx = 3, supp.type = "count", load2R = F)})
# GSE230792.cnt[1:5,1:5]
GSE189935.time = system.time({GSE189935.cnt = ParseGEO(acce = "GSE189935", down.supp = TRUE, supp.idx = 3, supp.type = "count", load2R = F)})
# GSE189935.cnt[1:5,1:5]
GSE239582.time = system.time({GSE239582.cnt = ParseGEO(acce = "GSE239582", down.supp = TRUE, supp.idx = 3, supp.type = "count", load2R = F)})
# GSE239582.cnt[1:5,1:5]
GSE225891.time = system.time({GSE225891.cnt = ParseGEO(acce = "GSE225891", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE225891.cnt[1:5,1:5]
GSE181709.time = system.time({GSE181709.cnt = ParseGEO(acce = "GSE181709", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE181709.cnt[1:5,1:5]
GSE244559.time = system.time({GSE244559.cnt = ParseGEO(acce = "GSE244559", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE244559.cnt[1:5,1:5]
GSE237838.time = system.time({GSE237838.cnt = ParseGEO(acce = "GSE237838", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE237838.cnt[1:5,1:5]
GSE192616.time = system.time({GSE192616.cnt = ParseGEO(acce = "GSE192616", down.supp = TRUE, supp.idx = 3, file.regex = "ReadsPerGene.out.tab", supp.type = "count", load2R = F)}) # mixture
# GSE192616.cnt[1:5,1:5]
GSE195839.time = system.time({GSE195839.cnt = ParseGEO(acce = "GSE195839", down.supp = TRUE, supp.idx = 1, supp.type = "count", load2R = F)})
# GSE195839.cnt[1:5,1:5]

save.image("GEO_counts.RData")
rm(list = ls())

# * 10x Genomics, DNBelab C4, SeekOne, MobiDrop ---------------
GSE271836.time = system.time({GSE271836.cnt = ParseGEO(acce = "GSE271836", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE292908.time = system.time({GSE292908.cnt = ParseGEO(acce = "GSE292908", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE275038.time = system.time({GSE275038.cnt = ParseGEO(acce = "GSE275038", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE305141.time = system.time({GSE305141.cnt = ParseGEO(acce = "GSE305141", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE301654.time = system.time({GSE301654.cnt = ParseGEO(acce = "GSE301654", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE285164.time = system.time({GSE285164.cnt = ParseGEO(acce = "GSE285164", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE264607.time = system.time({GSE264607.cnt = ParseGEO(acce = "GSE264607", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE262964.time = system.time({GSE262964.cnt = ParseGEO(acce = "GSE262964", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE237524.time = system.time({GSE237524.cnt = ParseGEO(acce = "GSE237524", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE311766.time = system.time({GSE311766.cnt = ParseGEO(acce = "GSE311766", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE306904.time = system.time({GSE306904.cnt = ParseGEO(acce = "GSE306904", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE192616.time = system.time({GSE192616.cnt = ParseGEO(acce = "GSE192616", down.supp = TRUE, supp.idx = 3, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE304316.time = system.time({GSE304316.cnt = ParseGEO(acce = "GSE304316", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE289417.time = system.time({GSE289417.cnt = ParseGEO(acce = "GSE289417", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE289082.time = system.time({GSE289082.cnt = ParseGEO(acce = "GSE289082", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE311149.time = system.time({GSE311149.cnt = ParseGEO(acce = "GSE311149", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE306678.time = system.time({GSE306678.cnt = ParseGEO(acce = "GSE306678", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE274058.time = system.time({GSE274058.cnt = ParseGEO(acce = "GSE274058", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE271622.time = system.time({GSE271622.cnt = ParseGEO(acce = "GSE271622", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE250411.time = system.time({GSE250411.cnt = ParseGEO(acce = "GSE250411", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE207991.time = system.time({GSE207991.cnt = ParseGEO(acce = "GSE207991", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE248507.time = system.time({GSE248507.cnt = ParseGEO(acce = "GSE248507", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE311407.time = system.time({GSE311407.cnt = ParseGEO(acce = "GSE311407", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE311825.time = system.time({GSE311825.cnt = ParseGEO(acce = "GSE311825", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE291575.time = system.time({GSE291575.cnt = ParseGEO(acce = "GSE291575", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE275680.time = system.time({GSE275680.cnt = ParseGEO(acce = "GSE275680", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE310060.time = system.time({GSE310060.cnt = ParseGEO(acce = "GSE310060", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE303003.time = system.time({GSE303003.cnt = ParseGEO(acce = "GSE303003", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE292493.time = system.time({GSE292493.cnt = ParseGEO(acce = "GSE292493", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE281009.time = system.time({GSE281009.cnt = ParseGEO(acce = "GSE281009", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE277137.time = system.time({GSE277137.cnt = ParseGEO(acce = "GSE277137", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE240460.time = system.time({GSE240460.cnt = ParseGEO(acce = "GSE240460", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE310574.time = system.time({GSE310574.cnt = ParseGEO(acce = "GSE310574", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE279475.time = system.time({GSE279475.cnt = ParseGEO(acce = "GSE279475", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE253859.time = system.time({GSE253859.cnt = ParseGEO(acce = "GSE253859", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE282768.time = system.time({GSE282768.cnt = ParseGEO(acce = "GSE282768", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE289110.time = system.time({GSE289110.cnt = ParseGEO(acce = "GSE289110", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE310290.time = system.time({GSE310290.cnt = ParseGEO(acce = "GSE310290", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE293930.time = system.time({GSE293930.cnt = ParseGEO(acce = "GSE293930", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE293928.time = system.time({GSE293928.cnt = ParseGEO(acce = "GSE293928", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE291111.time = system.time({GSE291111.cnt = ParseGEO(acce = "GSE291111", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE280910.time = system.time({GSE280910.cnt = ParseGEO(acce = "GSE280910", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE302089.time = system.time({GSE302089.cnt = ParseGEO(acce = "GSE302089", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE295144.time = system.time({GSE295144.cnt = ParseGEO(acce = "GSE295144", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE298041.time = system.time({GSE298041.cnt = ParseGEO(acce = "GSE298041", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE310026.time = system.time({GSE310026.cnt = ParseGEO(acce = "GSE310026", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE300217.time = system.time({GSE300217.cnt = ParseGEO(acce = "GSE300217", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE301305.time = system.time({GSE301305.cnt = ParseGEO(acce = "GSE301305", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE302421.time = system.time({GSE302421.cnt = ParseGEO(acce = "GSE302421", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE278892.time = system.time({GSE278892.cnt = ParseGEO(acce = "GSE278892", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE218499.time = system.time({GSE218499.cnt = ParseGEO(acce = "GSE218499", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE306115.time = system.time({GSE306115.cnt = ParseGEO(acce = "GSE306115", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE293929.time = system.time({GSE293929.cnt = ParseGEO(acce = "GSE293929", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE289049.time = system.time({GSE289049.cnt = ParseGEO(acce = "GSE289049", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE270070.time = system.time({GSE270070.cnt = ParseGEO(acce = "GSE270070", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE294062.time = system.time({GSE294062.cnt = ParseGEO(acce = "GSE294062", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE277570.time = system.time({GSE277570.cnt = ParseGEO(acce = "GSE277570", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE304467.time = system.time({GSE304467.cnt = ParseGEO(acce = "GSE304467", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE300477.time = system.time({GSE300477.cnt = ParseGEO(acce = "GSE300477", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE275148.time = system.time({GSE275148.cnt = ParseGEO(acce = "GSE275148", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE304771.time = system.time({GSE304771.cnt = ParseGEO(acce = "GSE304771", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE304846.time = system.time({GSE304846.cnt = ParseGEO(acce = "GSE304846", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE268811.time = system.time({GSE268811.cnt = ParseGEO(acce = "GSE268811", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE282929.time = system.time({GSE282929.cnt = ParseGEO(acce = "GSE282929", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE274844.time = system.time({GSE274844.cnt = ParseGEO(acce = "GSE274844", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE294399.time = system.time({GSE294399.cnt = ParseGEO(acce = "GSE294399", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})
GSE288941.time = system.time({GSE288941.cnt = ParseGEO(acce = "GSE288941", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", load2R = F, timeout = 36000)})

save.image("GEO_10x.RData")
rm(list = ls())

# Processed object -----------------------
# * rds, rds.gz -----------------
GSE285723.time = system.time({ParseGEOProcessed(acce = "GSE285723", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE311407.time = system.time({ParseGEOProcessed(acce = "GSE311407", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE271136.time = system.time({ParseGEOProcessed(acce = "GSE271136", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE270894.time = system.time({ParseGEOProcessed(acce = "GSE270894", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE306192.time = system.time({ParseGEOProcessed(acce = "GSE306192", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE280954.time = system.time({ParseGEOProcessed(acce = "GSE280954", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE299152.time = system.time({ParseGEOProcessed(acce = "GSE299152", timeout = 360000, supp.idx = 8, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE298041.time = system.time({ParseGEOProcessed(acce = "GSE298041", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE304675.time = system.time({ParseGEOProcessed(acce = "GSE304675", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE297041.time = system.time({ParseGEOProcessed(acce = "GSE297041", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE270486.time = system.time({ParseGEOProcessed(acce = "GSE270486", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE287517.time = system.time({ParseGEOProcessed(acce = "GSE287517", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE300712.time = system.time({ParseGEOProcessed(acce = "GSE300712", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE286410.time = system.time({ParseGEOProcessed(acce = "GSE286410", timeout = 360000, supp.idx = 3, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE287112.time = system.time({ParseGEOProcessed(acce = "GSE287112", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE272237.time = system.time({ParseGEOProcessed(acce = "GSE272237", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE309625.time = system.time({ParseGEOProcessed(acce = "GSE309625", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE308848.time = system.time({ParseGEOProcessed(acce = "GSE308848", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE301185.time = system.time({ParseGEOProcessed(acce = "GSE301185", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE299342.time = system.time({ParseGEOProcessed(acce = "GSE299342", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE309017.time = system.time({ParseGEOProcessed(acce = "GSE309017", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE307976.time = system.time({ParseGEOProcessed(acce = "GSE307976", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE274544.time = system.time({ParseGEOProcessed(acce = "GSE274544", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE288912.time = system.time({ParseGEOProcessed(acce = "GSE288912", timeout = 360000, supp.idx = 6, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE307539.time = system.time({ParseGEOProcessed(acce = "GSE307539", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})

# * RData, RData.gz ---------------
GSE311825.time = system.time({ParseGEOProcessed(acce = "GSE311825", timeout = 360000, supp.idx = 4, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE277678.time = system.time({ParseGEOProcessed(acce = "GSE277678", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE306745.time = system.time({ParseGEOProcessed(acce = "GSE306745", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE293701.time = system.time({ParseGEOProcessed(acce = "GSE293701", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE286398.time = system.time({ParseGEOProcessed(acce = "GSE286398", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE241573.time = system.time({ParseGEOProcessed(acce = "GSE241573", timeout = 360000, supp.idx = 5, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE249307.time = system.time({ParseGEOProcessed(acce = "GSE249307", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE244572.time = system.time({ParseGEOProcessed(acce = "GSE244572", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE294410.time = system.time({ParseGEOProcessed(acce = "GSE294410", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE255455.time = system.time({ParseGEOProcessed(acce = "GSE255455", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE255103.time = system.time({ParseGEOProcessed(acce = "GSE255103", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE264104.time = system.time({ParseGEOProcessed(acce = "GSE264104", timeout = 360000, supp.idx = 3, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE291342.time = system.time({ParseGEOProcessed(acce = "GSE291342", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE282783.time = system.time({ParseGEOProcessed(acce = "GSE282783", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})

# * h5ad, h5ad.gz ---------------
GSE311813.time = system.time({ParseGEOProcessed(acce = "GSE311813", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE311811.time = system.time({ParseGEOProcessed(acce = "GSE311811", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE309972.time = system.time({ParseGEOProcessed(acce = "GSE309972", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE295484.time = system.time({ParseGEOProcessed(acce = "GSE295484", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE308968.time = system.time({ParseGEOProcessed(acce = "GSE308968", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE289791.time = system.time({ParseGEOProcessed(acce = "GSE289791", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE288481.time = system.time({ParseGEOProcessed(acce = "GSE288481", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE250516.time = system.time({ParseGEOProcessed(acce = "GSE250516", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE305979.time = system.time({ParseGEOProcessed(acce = "GSE305979", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE285972.time = system.time({ParseGEOProcessed(acce = "GSE285972", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE290433.time = system.time({ParseGEOProcessed(acce = "GSE290433", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE290106.time = system.time({ParseGEOProcessed(acce = "GSE290106", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE288653.time = system.time({ParseGEOProcessed(acce = "GSE288653", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE262225.time = system.time({ParseGEOProcessed(acce = "GSE262225", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE277678.time = system.time({ParseGEOProcessed(acce = "GSE277678", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE297835.time = system.time({ParseGEOProcessed(acce = "GSE297835", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE307295.time = system.time({ParseGEOProcessed(acce = "GSE307295", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})

# * loom, loom.gz ---------------
GSE286325.time = system.time({ParseGEOProcessed(acce = "GSE286325", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE305371.time = system.time({ParseGEOProcessed(acce = "GSE305371", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE240922.time = system.time({ParseGEOProcessed(acce = "GSE240922", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE282790.time = system.time({ParseGEOProcessed(acce = "GSE282790", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE290839.time = system.time({ParseGEOProcessed(acce = "GSE290839", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE289737.time = system.time({ParseGEOProcessed(acce = "GSE289737", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE241403.time = system.time({ParseGEOProcessed(acce = "GSE241403", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE196796.time = system.time({ParseGEOProcessed(acce = "GSE196796", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE250597.time = system.time({ParseGEOProcessed(acce = "GSE250597", timeout = 360000, supp.idx = 4, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE254944.time = system.time({ParseGEOProcessed(acce = "GSE254944", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE247509.time = system.time({ParseGEOProcessed(acce = "GSE247509", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE256395.time = system.time({ParseGEOProcessed(acce = "GSE256395", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE246714.time = system.time({ParseGEOProcessed(acce = "GSE246714", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})

save.image("GEO_processed.RData")


