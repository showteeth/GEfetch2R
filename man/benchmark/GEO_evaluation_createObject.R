library(tidyverse)
library(GEfetch2R)

# Count matrix -------------
# * Smart-Seq2 ---------------
GSE94820.time = system.time({GSE94820.seu = ParseGEO(acce = "GSE94820", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE94820.seu
GSE214129.time = system.time({GSE214129.seu = ParseGEO(acce = "GSE214129", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE214129.seu
GSE292261.time = system.time({GSE292261.seu = ParseGEO(acce = "GSE292261", down.supp = TRUE, supp.idx = 2, supp.type = "count", timeout = 36000)})
# GSE292261.seu
GSE224915.time = system.time({GSE224915.seu = ParseGEO(acce = "GSE224915", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE224915.seu
GSE185076.time = system.time({GSE185076.seu = ParseGEO(acce = "GSE185076", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE185076.seu
GSE241004.time = system.time({GSE241004.seu = ParseGEO(acce = "GSE241004", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE241004.seu
GSE244538.time = system.time({GSE244538.seu = ParseGEO(acce = "GSE244538", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE244538.seu
GSE266889.time = system.time({GSE266889.seu = ParseGEO(acce = "GSE266889", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE266889.seu
GSE297431.time = system.time({GSE297431.seu = ParseGEO(acce = "GSE297431", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE297431.seu
GSE286046.time = system.time({GSE286046.seu = ParseGEO(acce = "GSE286046", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE286046.seu
GSE262789.time = system.time({GSE262789.seu = ParseGEO(acce = "GSE262789", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE262789.seu
GSE276332.time = system.time({GSE276332.seu = ParseGEO(acce = "GSE276332", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE276332.seu
GSE274944.time = system.time({GSE274944.seu = ParseGEO(acce = "GSE274944", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE274944.seu
GSE266915.time = system.time({GSE266915.seu = ParseGEO(acce = "GSE266915", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE266915.seu
GSE182219.time = system.time({GSE182219.seu = ParseGEO(acce = "GSE182219", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE182219.seu
GSE225056.time = system.time({GSE225056.seu = ParseGEO(acce = "GSE225056", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE225056.seu
GSE280276.time = system.time({GSE280276.seu = ParseGEO(acce = "GSE280276", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE280276.seu
GSE268201.time = system.time({GSE268201.seu = ParseGEO(acce = "GSE268201", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)}) # with comment char: '#
# GSE268201.seu
GSE252076.time = system.time({GSE252076.seu = ParseGEO(acce = "GSE252076", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE252076.seu
GSE284380.time = system.time({GSE284380.seu = ParseGEO(acce = "GSE284380", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE284380.seu
GSE245172.time = system.time({GSE245172.seu = ParseGEO(acce = "GSE245172", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE245172.seu
GSE215960.time = system.time({GSE215960.seu = ParseGEO(acce = "GSE215960", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE215960.seu
GSE261729.time = system.time({GSE261729.seu = ParseGEO(acce = "GSE261729", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE261729.seu
GSE266449.time = system.time({GSE266449.seu = ParseGEO(acce = "GSE266449", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE266449.seu
GSE264105.time = system.time({GSE264105.seu = ParseGEO(acce = "GSE264105", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)}) # anomalous data
# GSE264105.seu
GSE249746.time = system.time({GSE249746.seu = ParseGEO(acce = "GSE249746", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE249746.seu
GSE268974.time = system.time({GSE268974.seu = ParseGEO(acce = "GSE268974", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE268974.seu
GSE241165.time = system.time({GSE241165.seu = ParseGEO(acce = "GSE241165", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE241165.seu
GSE271100.time = system.time({GSE271100.seu = ParseGEO(acce = "GSE271100", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE271100.seu
GSE264553.time = system.time({GSE264553.seu = ParseGEO(acce = "GSE264553", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE264553.seu
GSE267774.time = system.time({GSE267774.seu = ParseGEO(acce = "GSE267774", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE267774.seu
GSE261258.time = system.time({GSE261258.seu = ParseGEO(acce = "GSE261258", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE261258.seu
GSE268038.time = system.time({GSE268038.seu = ParseGEO(acce = "GSE268038", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE268038.seu
GSE234572.time = system.time({GSE234572.seu = ParseGEO(acce = "GSE234572", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE234572.seu
GSE265909.time = system.time({GSE265909.seu = ParseGEO(acce = "GSE265909", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE265909.seu
GSE237867.time = system.time({GSE237867.seu = ParseGEO(acce = "GSE237867", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE237867.seu
GSE256466.time = system.time({GSE256466.seu = ParseGEO(acce = "GSE256466", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE256466.seu
GSE211550.time = system.time({GSE211550.seu = ParseGEO(acce = "GSE211550", down.supp = TRUE, supp.idx = 2, supp.type = "count", timeout = 36000)})
# GSE211550.seu
GSE207479.time = system.time({GSE207479.seu = ParseGEO(acce = "GSE207479", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE207479.seu
GSE241089.time = system.time({GSE241089.seu = ParseGEO(acce = "GSE241089", down.supp = TRUE, supp.idx = 5, supp.type = "count", timeout = 36000)})
# GSE241089.seu
GSE230792.time = system.time({GSE230792.seu = ParseGEO(acce = "GSE230792", down.supp = TRUE, supp.idx = 3, supp.type = "count", timeout = 36000)})
# GSE230792.seu
GSE189935.time = system.time({GSE189935.seu = ParseGEO(acce = "GSE189935", down.supp = TRUE, supp.idx = 3, supp.type = "count", timeout = 36000)})
# GSE189935.seu
GSE239582.time = system.time({GSE239582.seu = ParseGEO(acce = "GSE239582", down.supp = TRUE, supp.idx = 3, supp.type = "count", timeout = 36000)})
# GSE239582.seu
GSE225891.time = system.time({GSE225891.seu = ParseGEO(acce = "GSE225891", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE225891.seu
GSE181709.time = system.time({GSE181709.seu = ParseGEO(acce = "GSE181709", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE181709.seu
GSE244559.time = system.time({GSE244559.seu = ParseGEO(acce = "GSE244559", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE244559.seu
GSE237838.time = system.time({GSE237838.seu = ParseGEO(acce = "GSE237838", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE237838.seu
GSE192616.time = system.time({GSE192616.seu = ParseGEO(acce = "GSE192616", down.supp = TRUE, supp.idx = 3, file.regex = "ReadsPerGene.out.tab", supp.type = "count", timeout = 36000)}) # mixture
# GSE192616.seu
GSE195839.time = system.time({GSE195839.seu = ParseGEO(acce = "GSE195839", down.supp = TRUE, supp.idx = 1, supp.type = "count", timeout = 36000)})
# GSE195839.seu

save.image("GEO_counts_object.RData")
rm(list = ls())

# * 10x Genomics, SeekOne, MobiDrop ---------------
GSE271836.time = system.time({GSE271836.seu = ParseGEO(acce = "GSE271836", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE271836.seu
GSE292908.time = system.time({GSE292908.seu = ParseGEO(acce = "GSE292908", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE292908.seu
GSE275038.time = system.time({GSE275038.seu = ParseGEO(acce = "GSE275038", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE275038.seu
GSE305141.time = system.time({GSE305141.seu = ParseGEO(acce = "GSE305141", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE305141.seu
GSE301654.time = system.time({GSE301654.seu = ParseGEO(acce = "GSE301654", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE301654.seu
# GSE285164.time = system.time({GSE285164.seu = ParseGEO(acce = "GSE285164", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)}) # fail
# GSE285164.seu
GSE264607.time = system.time({GSE264607.seu = ParseGEO(acce = "GSE264607", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE264607.seu
GSE262964.time = system.time({GSE262964.seu = ParseGEO(acce = "GSE262964", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE262964.seu
GSE237524.time = system.time({GSE237524.seu = ParseGEO(acce = "GSE237524", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE237524.seu
GSE311766.time = system.time({GSE311766.seu = ParseGEO(acce = "GSE311766", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE311766.seu
GSE306904.time = system.time({GSE306904.seu = ParseGEO(acce = "GSE306904", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE306904.seu
GSE192616.time = system.time({GSE192616.seu = ParseGEO(acce = "GSE192616", down.supp = TRUE, supp.idx = 3, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE192616.seu
GSE304316.time = system.time({GSE304316.seu = ParseGEO(acce = "GSE304316", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE304316.seu
GSE289417.time = system.time({GSE289417.seu = ParseGEO(acce = "GSE289417", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE289417.seu
GSE289082.time = system.time({GSE289082.seu = ParseGEO(acce = "GSE289082", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE289082.seu
GSE311149.time = system.time({GSE311149.seu = ParseGEO(acce = "GSE311149", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE311149.seu
GSE306678.time = system.time({GSE306678.seu = ParseGEO(acce = "GSE306678", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE306678.seu
GSE274058.time = system.time({GSE274058.seu = ParseGEO(acce = "GSE274058", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE274058.seu
GSE271622.time = system.time({GSE271622.seu = ParseGEO(acce = "GSE271622", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE271622.seu
GSE250411.time = system.time({GSE250411.seu = ParseGEO(acce = "GSE250411", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE250411.seu
GSE207991.time = system.time({GSE207991.seu = ParseGEO(acce = "GSE207991", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE207991.seu
GSE248507.time = system.time({GSE248507.seu = ParseGEO(acce = "GSE248507", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE248507.seu
GSE311407.time = system.time({GSE311407.seu = ParseGEO(acce = "GSE311407", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE311407.seu
GSE311825.time = system.time({GSE311825.seu = ParseGEO(acce = "GSE311825", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE311825.seu
GSE291575.time = system.time({GSE291575.seu = ParseGEO(acce = "GSE291575", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE291575.seu
GSE275680.time = system.time({GSE275680.seu = ParseGEO(acce = "GSE275680", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE275680.seu
GSE310060.time = system.time({GSE310060.seu = ParseGEO(acce = "GSE310060", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE310060.seu
GSE303003.time = system.time({GSE303003.seu = ParseGEO(acce = "GSE303003", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE303003.seu
GSE292493.time = system.time({GSE292493.seu = ParseGEO(acce = "GSE292493", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE292493.seu
GSE281009.time = system.time({GSE281009.seu = ParseGEO(acce = "GSE281009", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE281009.seu
GSE277137.time = system.time({GSE277137.seu = ParseGEO(acce = "GSE277137", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE277137.seu
GSE240460.time = system.time({GSE240460.seu = ParseGEO(acce = "GSE240460", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE240460.seu
GSE310574.time = system.time({GSE310574.seu = ParseGEO(acce = "GSE310574", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE310574.seu
GSE279475.time = system.time({GSE279475.seu = ParseGEO(acce = "GSE279475", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE279475.seu
GSE253859.time = system.time({GSE253859.seu = ParseGEO(acce = "GSE253859", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE253859.seu
GSE282768.time = system.time({GSE282768.seu = ParseGEO(acce = "GSE282768", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE282768.seu
GSE289110.time = system.time({GSE289110.seu = ParseGEO(acce = "GSE289110", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE289110.seu
GSE310290.time = system.time({GSE310290.seu = ParseGEO(acce = "GSE310290", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE310290.seu
GSE293930.time = system.time({GSE293930.seu = ParseGEO(acce = "GSE293930", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE293930.seu
GSE293928.time = system.time({GSE293928.seu = ParseGEO(acce = "GSE293928", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE293928.seu
GSE291111.time = system.time({GSE291111.seu = ParseGEO(acce = "GSE291111", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE291111.seu
GSE280910.time = system.time({GSE280910.seu = ParseGEO(acce = "GSE280910", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE280910.seu
GSE302089.time = system.time({GSE302089.seu = ParseGEO(acce = "GSE302089", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE302089.seu
GSE295144.time = system.time({GSE295144.seu = ParseGEO(acce = "GSE295144", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE295144.seu
GSE298041.time = system.time({GSE298041.seu = ParseGEO(acce = "GSE298041", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE298041.seu
GSE310026.time = system.time({GSE310026.seu = ParseGEO(acce = "GSE310026", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE310026.seu
GSE300217.time = system.time({GSE300217.seu = ParseGEO(acce = "GSE300217", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE300217.seu
GSE301305.time = system.time({GSE301305.seu = ParseGEO(acce = "GSE301305", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE301305.seu
GSE302421.time = system.time({GSE302421.seu = ParseGEO(acce = "GSE302421", down.supp = TRUE, supp.idx = 1, supp.type = "10x", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE302421.seu
GSE278892.time = system.time({GSE278892.seu = ParseGEO(acce = "GSE278892", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE278892.seu
GSE218499.time = system.time({GSE218499.seu = ParseGEO(acce = "GSE218499", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE218499.seu
GSE306115.time = system.time({GSE306115.seu = ParseGEO(acce = "GSE306115", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE306115.seu
GSE293929.time = system.time({GSE293929.seu = ParseGEO(acce = "GSE293929", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE293929.seu
GSE289049.time = system.time({GSE289049.seu = ParseGEO(acce = "GSE289049", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE289049.seu
GSE270070.time = system.time({GSE270070.seu = ParseGEO(acce = "GSE270070", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE270070.seu
GSE294062.time = system.time({GSE294062.seu = ParseGEO(acce = "GSE294062", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE294062.seu
GSE277570.time = system.time({GSE277570.seu = ParseGEO(acce = "GSE277570", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE277570.seu
GSE304467.time = system.time({GSE304467.seu = ParseGEO(acce = "GSE304467", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE304467.seu
GSE300477.time = system.time({GSE300477.seu = ParseGEO(acce = "GSE300477", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE300477.seu
GSE275148.time = system.time({GSE275148.seu = ParseGEO(acce = "GSE275148", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE275148.seu
GSE304771.time = system.time({GSE304771.seu = ParseGEO(acce = "GSE304771", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE304771.seu
GSE304846.time = system.time({GSE304846.seu = ParseGEO(acce = "GSE304846", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE304846.seu
GSE268811.time = system.time({GSE268811.seu = ParseGEO(acce = "GSE268811", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE268811.seu
GSE282929.time = system.time({GSE282929.seu = ParseGEO(acce = "GSE282929", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE282929.seu
GSE274844.time = system.time({GSE274844.seu = ParseGEO(acce = "GSE274844", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE274844.seu
GSE294399.time = system.time({GSE294399.seu = ParseGEO(acce = "GSE294399", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE294399.seu
GSE288941.time = system.time({GSE288941.seu = ParseGEO(acce = "GSE288941", down.supp = TRUE, supp.type = "10xSingle", out.folder = "/home/songyabing/tmp/scfetch/GEOBench", merge = FALSE, timeout = 36000)})
GSE288941.seu

save.image("GEO_10x_object.RData")
rm(list = ls())

# Processed object -----------------------
# * rds, rds.gz -----------------
GSE285723.time = system.time({GSE285723.seu = ParseGEOProcessed(acce = "GSE285723", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE285723.seu
GSE311407.time = system.time({GSE311407.seu = ParseGEOProcessed(acce = "GSE311407", merge = FALSE, timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE311407.seu
GSE271136.time = system.time({GSE271136.seu = ParseGEOProcessed(acce = "GSE271136", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE271136.seu
GSE270894.time = system.time({GSE270894.seu = ParseGEOProcessed(acce = "GSE270894", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE270894.seu
GSE306192.time = system.time({GSE306192.seu = ParseGEOProcessed(acce = "GSE306192", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE306192.seu
# GSE280954.time = system.time({GSE280954.seu = ParseGEOProcessed(acce = "GSE280954", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)}) # fail with Seurat v4, success with Seurat v5
# GSE280954.seu
GSE299152.time = system.time({GSE299152.seu = ParseGEOProcessed(acce = "GSE299152", merge = FALSE, timeout = 360000, supp.idx = 8, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE299152.seu
GSE298041.time = system.time({GSE298041.seu = ParseGEOProcessed(acce = "GSE298041", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE298041.seu
GSE304675.time = system.time({GSE304675.seu = ParseGEOProcessed(acce = "GSE304675", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE304675.seu
GSE297041.time = system.time({GSE297041.seu = ParseGEOProcessed(acce = "GSE297041", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE297041.seu
GSE270486.time = system.time({GSE270486.seu = ParseGEOProcessed(acce = "GSE270486", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE270486.seu
GSE287517.time = system.time({GSE287517.seu = ParseGEOProcessed(acce = "GSE287517", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE287517.seu
GSE300712.time = system.time({GSE300712.seu = ParseGEOProcessed(acce = "GSE300712", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE300712.seu
GSE286410.time = system.time({GSE286410.seu = ParseGEOProcessed(acce = "GSE286410", merge = FALSE, timeout = 360000, supp.idx = 3, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE286410.seu
GSE287112.time = system.time({GSE287112.seu = ParseGEOProcessed(acce = "GSE287112", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE287112.seu
GSE272237.time = system.time({GSE272237.seu = ParseGEOProcessed(acce = "GSE272237", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE272237.seu
GSE309625.time = system.time({GSE309625.seu = ParseGEOProcessed(acce = "GSE309625", merge = FALSE, timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE309625.seu
GSE308848.time = system.time({GSE308848.seu = ParseGEOProcessed(acce = "GSE308848", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE308848.seu
GSE301185.time = system.time({GSE301185.seu = ParseGEOProcessed(acce = "GSE301185", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE301185.seu
GSE299342.time = system.time({GSE299342.seu = ParseGEOProcessed(acce = "GSE299342", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE299342.seu
# GSE309017.time = system.time({GSE309017.seu = ParseGEOProcessed(acce = "GSE309017", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)}) # fail with Seurat v4, success with Seurat v5
# GSE309017.seu
# GSE307976.time = system.time({GSE307976.seu = ParseGEOProcessed(acce = "GSE307976", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)}) # fail with Seurat v4, success with Seurat v5
# GSE307976.seu
GSE274544.time = system.time({GSE274544.seu = ParseGEOProcessed(acce = "GSE274544", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE274544.seu
GSE288912.time = system.time({GSE288912.seu = ParseGEOProcessed(acce = "GSE288912", merge = FALSE, timeout = 360000, supp.idx = 6, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE288912.seu
GSE307539.time = system.time({GSE307539.seu = ParseGEOProcessed(acce = "GSE307539", merge = FALSE, timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"), return.seu = T)})
GSE307539.seu

save.image("GEO_processed_object.RData")

# * Process RData, RData.gz files -----------------------------
GSE311825.time = system.time({ParseGEOProcessed(acce = "GSE311825", timeout = 360000, supp.idx = 4, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
# the count matrix can't be extracted correctly with Seurat v4, success with Seurat v5
GSE311825.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE311825/GSE311825_v3.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE277678.time = system.time({ParseGEOProcessed(acce = "GSE277678", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE277678.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE277678/GSE277678_AML_CR_Pre_Treatment_HMA-VEN_Refractory.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE306745.time = system.time({ParseGEOProcessed(acce = "GSE306745", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE306745.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE306745/GSE306745_Microglia2.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE293701.time = system.time({ParseGEOProcessed(acce = "GSE293701", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
# the count matrix can't be extracted correctly with Seurat v4, success with Seurat v5
GSE293701.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE293701/GSE293701_FINAL_B816-organoids.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE286398.time = system.time({ParseGEOProcessed(acce = "GSE286398", timeout = 360000, supp.idx = 2, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
# the count matrix can't be extracted correctly with Seurat v4, success with Seurat v5
GSE286398.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE286398/GSE286398_ctrl_ko_kobmt_ei_seur_objs.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE241573.time = system.time({ParseGEOProcessed(acce = "GSE241573", timeout = 360000, supp.idx = 5, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE241573.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE241573/GSE241573_Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE249307.time = system.time({ParseGEOProcessed(acce = "GSE249307", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE249307.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE249307/GSE249307_scRNA_seurat_data.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE244572.time = system.time({ParseGEOProcessed(acce = "GSE244572", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE244572.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE244572/GSE244572_RPE_CITESeq.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE294410.time = system.time({ParseGEOProcessed(acce = "GSE294410", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE294410.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE294410/GSE294410_rat_amygdala_obj.Rdata",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE255455.time = system.time({ParseGEOProcessed(acce = "GSE255455", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE255455.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE255455/GSM8072658_Sample1.Rdata",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE255103.time = system.time({ParseGEOProcessed(acce = "GSE255103", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE255103.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE255103/GSM8064141_Sample19.Rdata",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE264104.time = system.time({ParseGEOProcessed(acce = "GSE264104", timeout = 360000, supp.idx = 3, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE264104.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE264104/GSE264104_integrated_PBMC_cca_kanchor5_scClustViz.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE291342.time = system.time({ParseGEOProcessed(acce = "GSE291342", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
# the count matrix can't be extracted correctly with Seurat v4, success with Seurat v5
GSE291342.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE291342/GSE291342_snRNAseq_cellTypes_CLEAN.RData",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)
GSE282783.time = system.time({ParseGEOProcessed(acce = "GSE282783", timeout = 360000, supp.idx = 1, out.folder = "/home/songyabing/tmp/scfetch/GEOBench/Processed", file.ext = c("rdata", "rds", "h5ad", "loom"))})
GSE282783.list <- LoadRData(
  rdata = "/home/songyabing/tmp/scfetch/GEOBench/Processed/GSE282783/GSE282783_E16_FT_E17_hGFAP-Cre_mek12dcko.Rdata",
  accept.fmt = c("Seurat", "seurat", "SingleCellExperiment", "cell_data_set", "CellDataSet", "DESeqDataSet", "DGEList"),
  slot = "counts", return.obj = FALSE
)








