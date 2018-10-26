import compare

plotLog = True
parallel = True
compareFromPts = [32, 64,  128, 256, 512]
#compareToPts = [128,  256,  512]

#compareToPts = [128,  128,  128]
#compareToPts = [256, 256, 256, 256]
#compareToPts = [256, 256, 256]
compareToPts = [1024,  1024,  1024,  1024, 1024]
#compareToPts = [2048,  2048,  2048, 2048]

#compareFromPts = [32, 64, 128]

#compareFromPts = [32,  64,  128, 256,]
#compareToPts = [512,  512,  512,  512]

#where are the files for analysis located? end with /
dataDir = '/convection-in-sea-ice/test/bm1Diffusion/'
#outputDir = '/convection-in-sea-ice/test/bm1/oldModel/'

compare.convergenceStudy(compare.THETA, compareFromPts, compareToPts, compare.CALCULATED, compare.ANALYTIC, 
compare.L1ERR, plotLog, parallel, dataDir)
