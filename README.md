# Thesis-MS-data-processing
Overview of the MS data processing process. 
GUIDE - READ BEFORE PROCESSING
Content
This guide goes through the set up for processing from individual .mzML files to a a combined TIC sheet
It will go through the following aspects of processing:
1. Workspace
    What does the workspace need to contain?
2. Converting ‘.RAW’ to ‘.mzML’
    What will be assumed about your data
3. Metadata 
    The format it will need to be in
4.Processing settings
5. Sense checks

-----------------
1. Workspace
Do not work from the shared drive R-processing folder, create a new directory in your workspace. 
i.e. G:/Shared drives/Mass Spec _Metabolomics/James/NewExperimentProcessing
In this folder you will need: 
  A copy of the 'DI-ESI Processing Script' R-File 
  Your .mzML files in a folder called 'RAW' 
  An empty folder called 'SUM' 
  A 'Metadata.csv' file with your sample info on

Packages & R version
This script was written in R 4.10, newer versions may have changes that aren't compatible, so if there is a recurring error consider reverting back to this version. 
Packages will take some time to install, but you should only have to do it once (they will be saved to your library). During the installation process you may be asked if you want to install updates in the console, this may alter the functions in later versions so sticking with ‘n’ is a safe bet. Once installed you can just library them rather than having to reinstall
To this end the script has the installation code greyed out. You can remove the ‘#’ from the front and the code can be run if needed. 
Otherwise run the > lapply(c(“MSnbase”..... Code and it will library the required packages, and print ‘TRUE’ in the console if successful
If FALSE appears, reinstall the corresponding package
Firstly will need some mass spec data processing specific packages from a repository known as BiocManager (version 3.10 - later versions may work but this is what it was written with) 
  MSnbase
  MALDIquant
  MALDIquantForeign
  Dplyr (If you have a problem with the TechReps portion of the code it is likely because you haven't got dplyr installed and in the library)

---------------------------------
2. Converting ‘.RAW’ to ‘.mzML’
Prior to processing you need to convert your '.raw' files output by the mass spec into '.mzML' files
    This can be done using MSConvert (Proteowizard) 
    Here you can also cut out prewash/wash scans from your files 
        i.e. if your sample peak in the chromatogram is between scans 20 and 40, you can select only these scans
    Use the subset filter to do this, don't forget to press 'ADD'
      Check your .mzML files after they've been converted using SeeMS, look for:
        Should only have your selected scans in (although the entire chromatogram image will remain) 
        Should see the mass spectra seen in mass lynx on the actual MS computer
        Should all be similar size files (if there are some particular small/big ones double check them as it occasionally errors)
What will be assumed about your data
    The script will assume that every .mzML file in the raw folder wants processing, and that it is all in the metadata. 
    If it is not in the metadata it will stop running after the scan-summing step as the list of files will not match up.
    It will assume that all the scans in the .mzML files are wanted and will contribute to the end value
    If you only want a select scan range, specify this in the conversion step. 

---------------------------------------
3. Metadata
You will need a Metadata.csv file with all your sample information in it. 
It will need the first three columns to be:
  File 		[This is the file name the data is saved under] 
  SampleID	[This is the samples unique identifier] 
  Techrep 	[If only one taken still call it ‘1’]
  After this you can add any extra information you want and give it whatever header you want [If you do add new columns don't leave any file rows with blank cells]
Note: If you have acquired technical reps, the output will give you a table of masses in each rep and a separate table where technical reps are averaged.

----------------------------------------
4. Processing Settings
There are two sections of the processing that we can edit to fit our needs. Below will go through the basic considerations. 
NOTE: This script is written to give a very high resolution detection rate (good for noisy samples with low intensity) and a standard [low] resolution detection rate (Good for most normal experiments). 
Processing steps basic guide:
   A. Transformation
   B. Smoothing
   C. Baseline correction
   D. Normalisation
   E. Peak alignment
   F. Peak detection
   G.  Peak binning

A. Transformation
Reasoning for transformation: 
    To stabilise the variance in the data
    The detected intensity of an analyte can vary due to many factors in the detection process, meaning that there can be considerable variation between peaks. 
    By normalising the peaks we can make them more comparable, and give more accurate and precise measurements
Types of transformations available and which to use:
    Log transformation {“log”, “log2”,“log10”}
        Good for removing the impact of huge outliers
        Often over corrects, leaving too much noise
    Square root transformation {“sqrt”}
        Better choice for normally distributed data 
      This is the default in the script 

B. Smoothing
Reasoning for smoothing: 
    To reduce the noise in the data 
    This makes it easier to identify actual peaks
      HOWEVER overcorrection can smooth out small peaks or closely separated peaks.
      For example, if you're looking for two different metabolites at 191.0152 and 191.0538 these masses can be merged into one peak, giving one combined value
Types of smoothing available and which to use: 
  SavitzkyGolay {"SavitzkyGolay"} 
      Uses polynomial function to smooth the peaks
      Higher accuracy, because not only neighbouring peaks taken into account 
      Good for preserving peak shape
  MovingAverage {"MovingAverage"}
      Uses neighbouring peaks to smooth the peak
      Lover accuracy but less computationally expensive
      Does not preserve the shapes of the peaks
Other parameters
    Half window size (HWS)
      Determines how many data points are included in the smoothing
        Larger value = more smoothing
      If you’re trying to remove noise from your data a larger HWS 
      Default in script is set to 4

C. Baseline correction
Reasoning for baseline correction: 
    It removes the background noise from the spectra.
    Noise from the instrument/matrix/sample/etc will always be present
    Removing it makes the real peaks more obvious
    Artefacts due to noise are removed so you don't get false peaks
    The quantification/calibration of the peaks intensity is more accurate as it is being measured from the same baseline
    Peaks starting on a higher baseline may incorrectly be considered higher intensity, by making sure the baseline is corrected this doesn't happen. 
Methods of baseline correction available and which to use: 
   Statistics-sensitive Non-linear Iterative Peak-clipping algorithm {"SNIP"}
      Uses basic threshold to remove baseline noise
      Threshold slightly above the noise ratio of the data
      Little chance of introducing artefacts into the data
   TopHat {"TopHat"} 
      Uses moving average to remove baseline noise
      The neighbouring data has most influence on whether a point is removed, meaning it can falsely add artefacts if there is a constant signal 
      More effective at removing noise but can introduce artefacts
  ConvexHull {"ConvexHull"} 
      Removes a convex hull, more effective than above options
      Good at flattening consistent baseline shifts
      But not good at removing baselines that aren't smooth
  Median {"median"}
      Removes the median of the data 
      More effective than convex hull on data that isn't smooth
      But still may introduce artefacts
Other parameters
  Iterations (when using SNIP) / halfWindowSize (when using TopHat)
      The higher the iterations/HWS the more points included in the smoothing process
      For high resolution data this value may need to be lowered 
      Default for the script is 250

D. Normalisation/calibration
Reasoning for normalisation: 
    Signal intensity can be fairly arbitrary, so we need to normalise it into a more usable format
    TIC is easier to work with and more comparable between samples where average intensity is not. 
Types of normalisation available and which to use: 
  Total Ion Current {"TIC"}
      This normalises the intensity data to the entire spectrum, so that each peak’s intensity is a percentage of the overall signal 
  Probabilistic Quotient Normalization {"PQN"} 
      This divides the value for each peak by the signal of its precursor ion (the ion responsible for the peak) 
        The problem is it is majorly computationally expensive
        It is more accurate though as it will prevent over/under representation of certain compounds. 
  Median normalisation {"median"}
      Uses the median of the data to normalise to. 
      Less accurate but better at dealing with outliers. 

E. Peak alignment
Reasoning for peak alignment: 
      The detection of peaks to such fine masses means that there will be discrepancies between the m/z values for different samples. 
      These discrepancies between samples would mean inaccurate calling of a peak's intensity, when they are the same peak
      Essentially it makes sure you are comparing the same peaks. 
Parameters which to use: 
    Warping method = “Lowess” 
      Works well on noisy data and data with outliers 
      Based on local density of points
      There is the option for quadratic, linear and cubic but Lowess works best
  HalfWindowSize = 5 
      Similar to before it determines how many points are included in the ‘peak’ 
      If you want higher resolution (aka more separation) you can decrease this number 
  Signal to Noise ratio (SNR) = 3
      This separates peaks from noise	
      A local maximum (data point) has to be higher than this ratio to be considered a peak
      Lowering this value will lower the threshold for what is called a peak, meaning low intensity data is picked up
  Tolerance = 1e-4
      This is the local threshold for what is an identical peak. 
      Default is set to 0.0001 (4 decimal places)

F. Peak detection
Reasoning for peak detection: 
    This is where actual peaks are selected instead of signal curves
    We want the intensity value of a single m/z value not the entire peak shape. 
    Here we can significantly impact the resolution of our output by changing how the peaks are detected. 
    This is where the high resolution or low resolution settings come in. 
Types of peak detection available and which to use: 
  Mean absolute deviation {“MAD”} 
      It uses the median data points to create a threshold, and then anything above this threshold is considered a peak and not noise. 
      Using the median of the data makes this test more robust as outliers don’t affect the outcome.  
  Friedman's Super Smoother {“SuperSmoother”}
      Uses a moving average to work out what a peak is and what is the base line 
      Less robust as outliers can impact the average
Parameters to consider:
   HalfWindowSize = [Low Resolution] 10, [High resolution] 3
      Same as before this is how many neighbouring points are included in the assessment of a peak.
      Decreasing the value makes the peak detection  more specific. 
  Signal to noise ratio (SNR) = [Low Resolution] 4, [High resolution] 2
      Lower value will call more peaks, higher the opposite 
  refineMZ = “descendPeak”
      This defines how a peak is identified
      descendPeak sees the area of the peak as from the tip down until the measured signal increases again (aka the triangle of the peak) 
      Closely positioned peaks that have been over-smoothed will be merged into one data point as there will be no upward angle of the second peak
      This can be an issue worth checking 
  Signal percentage = 50 
      This determines how much of the peak area contributed to the m/z value calculated
      The default has it set to 50%, so half the peak area determines the value

G. Peak binning
Reasoning for peak binning: 
    This is for aligning all the m/z values so that there isn't discrepancies across the entire sample sets
    It is unifying the mz value across all samples
Types of peak binning available and which to use: 
  Strict {“strict”}
      Creates bins that cannot contain 2 or more peaks in the same sample
      Provides more distinct peaks
  Relaxed {“relaxed”}
      Allows multiple peaks in the same bin 
      Can capture broad peaks better
   Binning tolerance =  [Low Resolution] 0.05, [High resolution] 0.005
      Determines how big the bin will be
      Smaller value = more bins (and with strict, more possible peaks) 

-------------------------
5. Sense check
At the bottom of the script there is a list of plots that can be printed to check if the processing is doing what you expect it to. 
You will need to input a mass range of an expected metabolite in your sample. 
I chose 190.9-191.1 because I was expecting 2 peaks in my data in this range, and both should be present in all samples. 
Print the plots one by one and see how they change. 
You can edit your parameters at the top of the script to see how it changes things, but remember to reprocess the sections you are testing in the main body of the processing (you can just click through it all again) 
There is also a parameters output which will give you the m/z levels detected using the settings chosen. 
