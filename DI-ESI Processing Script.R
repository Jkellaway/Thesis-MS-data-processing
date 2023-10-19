    #                                                     -=-=-=-=-=-=-=-=-=-=-=-
    #                                     =-=-=-=-=-=-=-  DI-ESI MS data processing =-=-=-=-=-=-=-=-
    #                                                     -=-=-=-=-=-=-=-=-=-=-=-
    # Make sure you have your .mzML files in the right place and every  thing else set up (SEE: 'READ BEFORE PROCESSING.gdoc')

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
#                                                   ========={EDIT THIS SECTION}=========
    # Set working directory 
setwd("C:/Users/bop19jkp/Desktop/test")


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
#                                                   ---------{RUN THIS SECTION}--------
    # Load the packages you will need - (if previously installed just load them, the install infor is hashed out in case it is needed)
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("MSnbase")
#BiocManager::install("MALDIquant")
#BiocManager::install("MALDIquantForeign")
#BiocManager::install("dplyr")
lapply(c("MSnbase", "MALDIquant", "MALDIquantForeign", "dplyr"), require, character.only = TRUE)   

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
#                                                   ========={EDIT THIS SECTION}=========
    # Set data file name (this will be included in the output file name)
file <- "D_2 & D_3"

    # Set the parameters for processing       (Options below - explained in .gdoc)
Transformation_method <- "sqrt"               # [default = "sqrt"] "log", "log2", "log10"
Smoothing_method <- "SavitzkyGolay"           # [default = "SavitzkyGolay"] "MovingAverage" 
Smoothing_halfWindowSize <- 4                 # [default = 4] Increasing value = more smoothing, i.e. more noise removal
Baseline.Removal_method <- "SNIP"             # [default = "SNIP"] "TopHat", "ConvexHull", "median"
Baseline.Removal_iterations <- 250            # [default = 250] Higher value = more smoothing, for high res data you may want to lower 
Normalisation_method <- "TIC"                 # [default = "TIC"] "PQN", "Median"
Alignment_halfWindowSize <- 5                 # [default = 5] Increasing value = more separation between peaks (higher resolution, but more m/z discrepancy)
Alignment_SNR <- 3                            # [default = 3] Smaller value means lower intensity peaks considered actual peaks (can pick up noise as a result)
Alignment_tolerance <- 1e-4                   # [default = 1e-4] local threshold for identical peaks, smaller = stricter peaks

    # High and low res detection/binning settings - change as wanted for higher/lower resolution
HR.Detect_halfWindowSize <- 2                 # [default = 2] 
LR.Detect_halfWindowSize <- 10                # [default = 10] 

HR.Detect_SNR <- 2                            # [default = 2] 
LR.Detect_SNR <- 4                            # [default = 4] 

HR.Binning_tolerance <- 0.005                 # [default = 0.005] 
LR.Binning_tolerance <- 0.05                  # [default = 0.05] 

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#                                                   ---------{RUN THIS SECTION}--------
    # Sum scans
    # Create a function that will be used to the output the summed files 
createSumFilenames<-function(filename, out.folder="SUM"){                
  file.path(out.folder,sub("\\.mzML",".sum.mzML",basename(filename)))
  }
    # Create a function to do the summing
mergeScans<-function(filename, out.folder="SUM", tol.ppm=15){            
  filename.out<-createSumFilenames(filename, out.folder)
  if(file.exists(filename.out)){
    warning(paste("Not processing",filename,"as output file",filename.out,"already exists\n"))
  } else {
    profiledata<-readMSData(filename, centroided. = F, msLevel. = 1, mode = 'onDisk') 
    profiledata.merge<-MSnbase::combineSpectra(MSpectra(spectra(profiledata)), weighted=F,  timeDomain=T, ppm=tol.ppm, intensityFun=sum, unionPeaks=T)
    writeMSData(as(profiledata.merge, "MSnExp"), file=filename.out, outformat="mzml")
    message(paste0("Processed ",filename))
  }
}
infiles<-Sys.glob("RAW/*.mzML")                           # Create the list of files to do
invisible(sapply(infiles, mergeScans))                    # Apply the summing function (THIS WILL TAKE A LITTLE TIME TO COMPLETE) 

# A summed '.mzML' file should be in there for each sample - You can check these files again in See-MS

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    # Organise the metadata
MetaData <- read.csv("Metadata.csv", header = T)                          # Import the metadata
Metadata_Headers <- colnames(MetaData)                                    # Create a vector of column names
SumFiles<-as.data.frame(Sys.glob("SUM/*.sum.mzML"), stringsAsFactors = F) # Create a data frame with the summed file names listed 
colnames(SumFiles)[1] <-"Summed_Filename"
SumFiles <- cbind(SumFiles, MetaData)                                     # Merge the List of files and metadata
SumFiles$tech_id<-paste0(SumFiles$SampleID,".",SumFiles$TechRep)          # Create a unique ID for all samples 

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Processing steps
suppressWarnings(rawdata<-importMzMl(SumFiles$Summed_Filename, centroided = F, verbose=F))          # Import the data
spectra_transformed <- transformIntensity(rawdata,method=Transformation_method)                     # Transform data
spectra_smoothed <- smoothIntensity(spectra_transformed, method=Smoothing_method, halfWindowSize=Smoothing_halfWindowSize) # Smooth data
spectra_baselineCorrection <- removeBaseline(spectra_smoothed, method=Baseline.Removal_method, iterations=Baseline.Removal_iterations) # Baseline correction
spectra_calibrated <- calibrateIntensity(spectra_baselineCorrection, method=Normalisation_method)   # Normalise data by Total Ion count 
spectra_aligned <- alignSpectra(spectra_calibrated, halfWindowSize=Alignment_halfWindowSize, SNR=Alignment_SNR, tolerance=Alignment_tolerance, warpingMethod="lowess") # Perform a peak alignment 

    # Processing for individual files
HR.detect_peaks<-detectPeaks(spectra_aligned, method="MAD", halfWindowSize=HR.Detect_halfWindowSize, SNR=HR.Detect_SNR, refineMz="descendPeak", signalPercentage=50) 
HR.bin_peaks<-binPeaks(HR.detect_peaks, method="strict", tolerance = HR.Binning_tolerance)
HR.featureMatrix <- intensityMatrix(HR.bin_peaks, spectra_aligned) 
HR.intensities_out <- data.frame(SumFiles[Metadata_Headers], HR.featureMatrix, check.names = F, stringsAsFactors = F)
HR.intensities_out<-rbind(c("rt", rep(1, ncol(HR.intensities_out)-1)), HR.intensities_out) 
HR.intensities_out <- setNames(data.frame(t(HR.intensities_out[,-1])), HR.intensities_out[,1])
write.csv(HR.intensities_out, row.names = T, file = paste0(file,"-HighRes.Intensities.csv"))

LR.detect_peaks<-detectPeaks(spectra_aligned, method="MAD", halfWindowSize=LR.Detect_halfWindowSize, SNR=LR.Detect_SNR, refineMz="descendPeak", signalPercentage=50) 
LR.bin_peaks<-binPeaks(LR.detect_peaks, method="strict", tolerance = LR.Binning_tolerance)
LR.featureMatrix <- intensityMatrix(LR.bin_peaks, spectra_aligned) 
LR.intensities_out <- data.frame(SumFiles[Metadata_Headers], LR.featureMatrix, check.names = F, stringsAsFactors = F)
LR.intensities_out<-rbind(c("rt", rep(1, ncol(LR.intensities_out)-1)), LR.intensities_out) 
LR.intensities_out <- setNames(data.frame(t(LR.intensities_out[,-1])), LR.intensities_out[,1])
write.csv(LR.intensities_out, row.names = T, file = paste0(file,"-LowRes.Intensities.csv"))

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    # Averaging technical reps

    #gets mean for spectra based on sample id (which there are 3 tech reps of each)
avgSpectra <- averageMassSpectra(spectra_aligned, labels = SumFiles$SampleID, method = "mean") 
avgSumFiles <- filter(SumFiles, TechRep == 1) #lists first tech rep file of each biorep as sample reference

    # Processing for averaged files 
avg.HR.detect_peaks<-detectPeaks(avgSpectra, method="MAD", halfWindowSize=HR.Detect_halfWindowSize, SNR=HR.Detect_SNR, refineMz="descendPeak", signalPercentage=50) 
avg.HR.bin_peaks<-binPeaks(avg.HR.detect_peaks, method="strict", tolerance = HR.Binning_tolerance)
avg.HR.featureMatrix <- intensityMatrix(avg.HR.bin_peaks, avgSpectra)
avg.HR.intensities_out <- data.frame(avgSumFiles[Metadata_Headers], avg.HR.featureMatrix, check.names = F, stringsAsFactors = F)
avg.HR.intensities_out<-rbind(c("rt", rep(1, ncol(avg.HR.intensities_out)-1)), avg.HR.intensities_out) 
avg.HR.intensities_out <- setNames(data.frame(t(avg.HR.intensities_out[,-1])), avg.HR.intensities_out[,1]) 
write.csv(avg.HR.intensities_out, row.names = T, file = paste0(file,"-Avg_HighRes.Intensities.csv"))

avg.LR.detect_peaks<-detectPeaks(avgSpectra, method="MAD", halfWindowSize=LR.Detect_halfWindowSize, SNR=LR.Detect_SNR, refineMz="descendPeak", signalPercentage=50) 
avg.LR.bin_peaks<-binPeaks(avg.LR.detect_peaks, method="strict", tolerance = LR.Binning_tolerance)
avg.LR.featureMatrix <- intensityMatrix(avg.LR.bin_peaks, avgSpectra)
avg.LR.intensities_out <- data.frame(avgSumFiles[Metadata_Headers], avg.LR.featureMatrix, check.names = F, stringsAsFactors = F)
avg.LR.intensities_out<-rbind(c("rt", rep(1, ncol(avg.LR.intensities_out)-1)), avg.LR.intensities_out)
avg.LR.intensities_out <- setNames(data.frame(t(avg.LR.intensities_out[,-1])), avg.LR.intensities_out[,1]) 
write.csv(avg.LR.intensities_out, row.names = T, file = paste0(file,"-Avg_LowRes.Intensities.csv"))

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    # Create a text doc with the parameters used in it 
ParametersUsed <- data.frame(Process= c("Transformation method","Smoothing method","Smoothing HWS", "Baseline correction method", "Baseline correction iterations", 
                      "Normalisation method", "Alignment method", "Alignment HWS", "Alignment SNR", "Alignment Tolerance","Peak detection method", "High Res peak detection HWS",
                      "Low Res peak detection HWS","High Res peak detection SNR","Low Res peak detection SNR",
                      "High Res peak binning tolerance","Low Res peak binning tolerance", "-",
                      "# of HighRes m/z levels","# of LowRes m/z levels","# of Avg_HighRes m/z levels","# of Avg_LowRes m/z levels",
                      "# of HighRes samples","# of LowRes samples","# of Avg_HighRes samples","# of Avg_LowRes samples"), 
           Setting = c(Transformation_method, Smoothing_method, Smoothing_halfWindowSize, Baseline.Removal_method, Baseline.Removal_iterations,
                       Normalisation_method, "Lowess",Alignment_halfWindowSize, Alignment_SNR, Alignment_tolerance,"MAD", HR.Detect_halfWindowSize,
                       LR.Detect_halfWindowSize, HR.Detect_SNR, LR.Detect_SNR, HR.Binning_tolerance, LR.Binning_tolerance, "-",
                       (ncol(HR.featureMatrix)-1), (ncol(LR.featureMatrix)-1), (ncol(avg.HR.featureMatrix)-1),(ncol(avg.LR.featureMatrix)-1),
                       nrow(HR.featureMatrix), nrow(LR.featureMatrix), nrow(avg.HR.featureMatrix),nrow(avg.LR.featureMatrix)))
write.table(ParametersUsed, row.names = F, file = paste0(file, "- Parameters used in processing.txt"), sep = " | ")


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
#                                                   ========={EDIT THIS SECTION}=========
    # Sense checks
  
    # set x-axis limits around a peak you expect to see in your data 
expected_peak <- c(190.9,191.1)

    # Check the raw spectrum is what you expect
plot(rawdata[[1]], main = "Rawdata")
    
    # See how the transformation has lessens the impact of outliers to give better general signal  
plot(spectra_transformed[[1]], main = paste0("Transformed Data - ",Transformation_method))

    # see whether enough smoothing is occuring to remove the base noise
plot(spectra_transformed[[1]], main = paste0("Smoothed data - ",Smoothing_method,", HWS=", Smoothing_halfWindowSize), xlim =expected_peak) + 
  lines(spectra_smoothed[[1]], col = "blue")

    # See how the baseline is corrected
b <- estimateBaseline(spectra_smoothed[[1]], method = Baseline.Removal_method, iterations = Baseline.Removal_iterations)
plot(spectra_smoothed[[1]], main = "Baseline Correction needed") + lines(b, col = "red")
plot(spectra_baselineCorrection[[1]], main = paste0("Baseline Corrected - ", Baseline.Removal_method, ", iterations=", Baseline.Removal_iterations))

    # see how normalization changes intensity values
plot(spectra_calibrated[[1]], main = paste0("Normalisation by ",Normalisation_method))

    # See how aligning spectra is needed (might not change much)
plot(spectra_calibrated[[1]], main = "Not aligned" , xlim = expected_peak, type = "b") + 
  lines(spectra_calibrated[[5]], col = "blue", type = "b") + 
  lines(spectra_calibrated[[10]], col = "red", type = "b") + 
  lines(spectra_calibrated[[16]], col = "Gold", type = "b")
plot(spectra_aligned[[1]], main = paste0("Aligned; HWS=", Alignment_halfWindowSize, ", SNR=", Alignment_SNR, ", Tolerance=", Alignment_tolerance), xlim = expected_peak, type = "b") + 
  lines(spectra_aligned[[5]], col = "blue", type = "b") + 
  lines(spectra_aligned[[10]], col = "red", type = "b") + 
  lines(spectra_aligned[[16]], col = "Gold", type = "b")

    # Show how detecting peaks turns the curve into a point (high and low res)
plot(HR.detect_peaks[[1]], main = paste0("HR.Peak detection; HWS=", HR.Detect_halfWindowSize, ", SNR=",HR.Detect_SNR), xlim = expected_peak) +
  lines(spectra_aligned[[1]], col = "black") +
  lines(spectra_aligned[[5]], col = "blue") + 
  lines(spectra_aligned[[10]], col = "red") + 
  lines(spectra_aligned[[16]], col = "Gold") +  
  lines(HR.detect_peaks[[1]], col = "black") +
  lines(HR.detect_peaks[[5]], col = "blue") + 
  lines(HR.detect_peaks[[10]], col = "red") + 
  lines(HR.detect_peaks[[16]], col = "Gold")

plot(LR.detect_peaks[[1]], main = paste0("LR.Peak detection; HWS=", LR.Detect_halfWindowSize, ", SNR=",LR.Detect_SNR), xlim = expected_peak) +
  lines(spectra_aligned[[1]], col = "black") +
  lines(spectra_aligned[[5]], col = "blue") + 
  lines(spectra_aligned[[10]], col = "red") + 
  lines(spectra_aligned[[16]], col = "Gold") +  
  lines(LR.detect_peaks[[1]], col = "black") +
  lines(LR.detect_peaks[[5]], col = "blue") + 
  lines(LR.detect_peaks[[10]], col = "red") + 
  lines(LR.detect_peaks[[16]], col = "Gold")

    # see how peak binning effects the data (low and high res)
plot(HR.detect_peaks[[1]], main = "HR. peaks - no binning", xlim = expected_peak)+
  lines(HR.detect_peaks[[5]], col = "blue") + 
  lines(HR.detect_peaks[[10]], col = "red") + 
  lines(HR.detect_peaks[[16]], col = "Gold")
plot(HR.bin_peaks[[1]], main = paste0("HR. peaks - binned by ", HR.Binning_tolerance), xlim = expected_peak) +
  lines(HR.bin_peaks[[5]], col = "blue") + 
  lines(HR.bin_peaks[[10]], col = "red") + 
  lines(HR.bin_peaks[[16]], col = "Gold")

plot(LR.detect_peaks[[1]], main = "LR. peaks - no binning", xlim = expected_peak)+
  lines(LR.detect_peaks[[5]], col = "blue") + 
  lines(LR.detect_peaks[[10]], col = "red") + 
  lines(LR.detect_peaks[[16]], col = "Gold")
plot(LR.bin_peaks[[1]], main = paste0("LR. peaks - binned by ", LR.Binning_tolerance), xlim = expected_peak) +
  lines(LR.bin_peaks[[5]], col = "blue") + 
  lines(LR.bin_peaks[[10]], col = "red") + 
  lines(LR.bin_peaks[[16]], col = "Gold")


