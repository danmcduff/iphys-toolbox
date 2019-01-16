# iPhys: An Open  Non-Contact  Imaging-Based  Physiological  Measurement Toolbox
In the past few years a lot of attention has been given to methods for remotely measuring physiological signals using low-cost cameras.  Imaging PPG (iPPG) focuses on the measurement of volumetric changes in blood flow at distance from the body using imaging devices to capture changes in transmitted or reflected light. Imaging ballistocardiography (iBCG) typically leverages optical flow estimation to track the vertical motion of the head or body from a video sequence. Both iPPG and iBCG methods can be used to recover human vital signals.

This toolbox contains MATLAB implementations of a number of algorithms for non-contact physiological measurement. This will enable researchers to present results on their datasets using standard public implementations of the baseline methods with all parameters known. The toolbox includes implementations of many of the most commonly used baseline methods for imaging photplethysmography (iPPG) and image ballistocardiography (iBCG).

## Citation: ##

If you find this toolbox helpful and make use of it in your work please cite:

### iPhys: An Open Non-Contact Imaging-Based Physiological Measurement Toolbox ###

Daniel McDuff

ArXiv

Link: http://arxiv.org/abs/1901.04366
Cite as:	arXiv:1901.04366 [cs.CV]

We welcome suggestions or contact regarding the toolbox.  Please contact: damcduff@microsoft.com

## Background: ## 

Imaging photoplethysmography (iPPG) has developed as a method for capturing the BVP signal remotely using digital cameras and ambient light. Almost any digital camera (i.e., webcam or cellphone camera) is sufficiently sensitive to capture the pulse signal when the subject is close to the device.  Poh et al. (2010) proposed the use of Independent Component Analysis (ICA) to recover the pulse signal and enable a fully automated iPPG framework.

![Alt text](imgs/Imaging_PPG.png?raw=true "Imaging PPG pipeline.")

Methods inspired by optical models of the skin have helped advance the state-of-the art. The CHROM (2013) method uses a linear combination of the chrominance signals and makes the assumption of a standardized skin color profile to white-balance the video frames. 
The plane orthogonal to the skin (POS [2017]) algorithm assumes the presence of a pulsatile color space signal and posits that this will be orthogonal to the skin color space.

## Licensing: ##

The code is licensed with the MIT License (https://opensource.org/licenses/MIT) and the RAIL AI Licenses (https://www.licenses.ai/) that allows use of the software including rights to use, copy or modify it, with the exception of harmful applications as specified in: https://www.licenses.ai/open-source-license.

## Methods: ## 

The toolbox contains implementations of the several of the most commonly used methods for video-based physiological measurement. Below we give citations and links to the publications on which these implementations are based:

### GREEN CHANNEL - Verkruysse, Svaasand, Nelson (2008) ###

https://www.osapublishing.org/oe/viewmedia.cfm?uri=oe-16-26-21434&seq=0

Citation: Verkruysse, W., Svaasand, L. O., & Nelson, J. S. (2008). Remote plethysmographic imaging using ambient light. Optics express, 16(26), 21434-21445. DOI: 10.1364/OE.16.021434

### ICA - Poh, McDuff, Picard (2010) ###

https://www.osapublishing.org/oe/viewmedia.cfm?uri=oe-18-10-10762&seq=0

Citation: Poh, M. Z., McDuff, D. J., & Picard, R. W. (2010). Non-contact, automated cardiac pulse measurements using video imaging and blind source separation. Optics express, 18(10), 10762-10774. DOI: 10.1364/OE.18.010762

### CHROM - De Haan & Jeanne (2013) ###

http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.726.6643&rep=rep1&type=pdf

Citation: De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886. DOI: 10.1109/TBME.2013.2266196

### POS - Wang, den Brinker, Stuijk & de Haan (2017) ###

https://pure.tue.nl/ws/files/31563684/TBME_00467_2016_R1_preprint.pdf

Citation: Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote PPG. IEEE Transactions on Biomedical Engineering, 64(7), 1479-1491. DOI: 10.1109/TBME.2016.2609282

### BCG - Balakrishnan, Durand & Guttag (2013) ###

http://openaccess.thecvf.com/content_cvpr_2013/papers/Balakrishnan_Detecting_Pulse_from_2013_CVPR_paper.pdf

Citation: Balakrishnan, G., Durand, F., & Guttag, J. (2013). Detecting pulse from head motions in video. In Computer Vision and Pattern Recognition (CVPR), 2013 IEEE Conference on (pp. 3430-3437). IEEE. DOI: 10.1109/CVPR.2013.440


We would happily include additional implementations of iPPG or iBCG algorithms in this toolbox. If you would like to contribute an implementations of a new method please submit a pull request or email: damcduff@microsoft.com

## Usage: ##

The scripts have been created to allow easy processing of any video dataset.  We have also included example data in the "test_data" folder.

Each method takes the following inputs:

*       VideoFile         =       Video file path.

*       FS                =       Video framerate (fps).

*       StartTime         =       Timepoint at which to start process (default = 0 seconds).

*       Duration          =       Duration of the time window to process (default = 60 seconds).

*       ECGFile           =       File path to corresponding ECG data (.mat) file containing: 1) The waveform - ECGData.data, 2) The ECG sampling rate - ECGData.fs, 3) The ECG peak locations (in samples) - ECGData.peaks.

*       PPGFile           =       File path to corresponding PPG data (.mat) file containing: 1) The waveform - PPGData.data, 2) The PPG sampling rate - PPGData.fs, 3) The PPG peak locations (in samples) - PPGData.peaks.

*       PlotTF            =       Boolean to turn plotting results on or off.

And produces the following output:

*        BVP              =    The predicted blood volume pulse signal.

*        PR               =    Estimated Pulse Rate (PR) from BVP timeseries using peak in periodogram.

*        HR_ECG           =    Gold standard reference Heart Rate (HR) measured from the ECG timeseries.

*        PR_PPG           =    Pulse Rate measured from the PPG timeseries for the window.

*        SNR              =    Blood Volume Pulse Signal-to-Noise Ratio calculated based on the ECG HR frequency.


### Example: ###
DataDirectory           = [cd '\test_data\'];

VideoFile               = [DataDirectory 'video_example.mp4'];

FS                      = 30;

StartTime               = 0;

Duration                = 60;

ECGFile                 = [DataDirectory 'ECGData.mat'];

PPGFile                 = [DataDirectory 'PPGData.mat'];

PlotTF                  = false;


[BVP, PR, HR_ECG, PR_PPG, SNR] = ICA_POH(VideoFile, FS, StartTime, Duration, ECGFile, PPGFile, PlotTF)


### Data ###

Due to the size of the video example, [please download it here](https://drive.google.com/open?id=1oD4VbBD9ColSlbiIMEgxbvQ7LnXHPy1_) and add it to the "test_data" folder. 

The ECGData.mat file has three fields: 1) The waveform - ECGData.data, 2) The ECG sampling rate - ECGData.fs, 3) The ECG peak locations (in samples) - ECGData.peaks.

The PPGData.mat file has three fields: 1) The waveform - PPGData.data, 2) The PPG sampling rate - PPGData.fs, 3) The PPG peak locations (in samples) - PPGData.peaks.
