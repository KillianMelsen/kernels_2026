set.seed(1997)
M <- as.data.frame(readxl::read_xlsx("Briwecs/raw_data/135K.xlsx"))
rownames(M) <- M$probeset_id
M <- t(as.matrix(M[, -1]))
M[which(M == "---")] <- NA
M2 <- ASRgenomics::snp.recode(M = M, na.string = "NA")
M2 <- M2$Mrecode
saveRDS(M2, "Briwecs/raw_data/135K_recoded.rds")
gdata <- statgenGWAS::createGData(geno = M2)
gdata <- statgenGWAS::codeMarkers(gdata, MAF = 0.05)
K <- statgenGWAS::kinship(gdata$markers, method = "IBS")
saveRDS(K, "Briwecs/raw_data/K.rds")
