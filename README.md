# iPhys: An Open  Non-Contact  Imaging-Based  Physiological  Measurement
In the past few years a lot of attention has been given to methods for remotely measuring PPG using low-cost cameras.  Imaging PPG (iPPG) focuses on the measurement of volumetric changes in blood flow at distance from the body using imaging devices to capture changes in transmitted or relected light. 

This toolbox contains MATLAB implementations of a number of algorithms for iPPG analysis.  The toolbox includes implementations of many of the most commonly used baseline methods and was created to accompany the largest public iPPG dataset.

If you find this toolbox helpful and make use of it in your work please cite:

An Open Imaging-based Physiological Measurement Toolbox

Daniel McDuff and Ethan Blackford

## Background: ## 

Imaging photoplethysmography (iPPG) has developed as a method for capturing the BVP signal remotely using digital cameras and ambient light. Almost any digital camera (i.e., webcam or cellphone camera) is sufficiently sensitive to capture the pulse signal when the subject is close to the device.  Poh et al. (2010) proposed the use of Independent Component Analysis (ICA) to recover the pulse signal and enable a fully automated iPPG framework

![Alt text](imgs/Imaging_PPG.png?raw=true "Imaging PPG pipeline.")

Methods inspired by optical models of the skin have helped advance the state-of-the art. The CHROM (2013) method uses a linear combination of the chrominance signals and makes the assumption of a standardized skin color profile to white-balance the video frames. 
The plane orthogonal to the skin (POS [2017]) algorithm assumes the presence of a pulsatile color space signal and posits that this will be orthogonal to the skin color space.

## Methods: ## 

The toolbox contains implementations of the several of the most commonly used methods for video-based physiological measurement. Below we give citations and links to the publications on which these implementations are based:

### GREEN CHANNEL - Verkruysse, Svaasand, Nelson (2008) ###

https://www.osapublishing.org/DirectPDFAccess/9A6ECF3A-E3F8-74A3-6FAD721D7E164888_175396/oe-16-26-21434.pdf?da=1&id=175396&seq=0&mobile=no

Verkruysse, W., Svaasand, L. O., & Nelson, J. S. (2008). Remote plethysmographic imaging using ambient light. Optics express, 16(26), 21434-21445.

### ICA - Poh, McDuff, Picard (2010) ###

https://www.osapublishing.org/viewmedia.cfm?uri=oe-18-10-10762&seq=0&origin=search

Citation: Poh, M. Z., McDuff, D. J., & Picard, R. W. (2010). Non-contact, automated cardiac pulse measurements using video imaging and blind source separation. Optics express, 18(10), 10762-10774.

### CHROM - De Haan & Jeanne (2013) ###

http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.726.6643&rep=rep1&type=pdf


Citation: De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886.

### POS - Wang, den Brinker, Stuijk & de Haan (2017) ###

https://pure.tue.nl/ws/files/31563684/TBME_00467_2016_R1_preprint.pdf

Citation: Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG. IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491.

### BCG - Balakrishnan, Durand & Guttag (2013) ###

http://openaccess.thecvf.com/content_cvpr_2013/papers/Balakrishnan_Detecting_Pulse_from_2013_CVPR_paper.pdf

Citation: Balakrishnan, G., Durand, F., & Guttag, J. (2013). Detecting pulse from head motions in video. In Computer Vision and Pattern Recognition (CVPR), 2013 IEEE Conference on (pp. 3430-3437). IEEE.


We would happily include additional implementations of iPPG or iBCG algorithms in this toolbox. If you would like to contribute an implemntations of a new method please submit a pull request or email: damcduff@microsoft.com

## Usage: ##

The iPPG-toolbox accompanies the AFRL iPPG dataset (described below).  The scripts have been created to allow easy processing of the dataset but also to generalize to other iPPG data.

Each method takes the following inputs:

*       VideoFile               =       Video filename.

*       FS                      =       Framerate of Video;

*       StartTime               =       Timepoint at which to start process (default = 15);

*       Duration                =       Duration of the time window to process (default = 30 seconds).

*       ECGData                 =       Corresponding ECGData data file.

*       PPGData                 =       Corresponding PPGData data file.

*       PlotTF                  =       Boolean to turn plotting results on or off.


And produces the following output:

*        BVP                     =    The predicted blood volume pulse signal.

*        PR                      =    Estimated Pulse Rate (PR) from BVP timeseries using peak in periodogram.

*        HR_ECG                  =    Gold startdard reference Heart Rate (HR) measured from the ECG timeseries.

*        PR_PPG                  =     Pulse Rate measured from the PPG timeseries for the window.

*        SNR                     =    Blood Volume Pulse Signal-to-Noise Ratio calculated based on the ECG HR frequency.


### Example: ###

VideoFile = test_data/example_video.mp4;

FS = 30;

StartTime = 15;

Duration = 300;

ECGData = test_data/ECGData;

PPGData = test_data/PPGData;

PlotTF = False;

[BVP, PR, HR_ECG, PR_PPG, SNR] = ICA_POH(VideoFile, FS, StartTime, Duration, ECGData, PPGData, PlotTF)


### Data ###

Due to the size of the video example, [please download it here](https://drive.google.com/open?id=1oD4VbBD9ColSlbiIMEgxbvQ7LnXHPy1_) and add it to the "test_data" folder. 
